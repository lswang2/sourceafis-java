// Part of SourceAFIS: https://sourceafis.machinezoo.com
package com.machinezoo.sourceafis;

import static java.util.stream.Collectors.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import java.util.stream.*;
import javax.imageio.*;
import com.google.gson.*;
import com.machinezoo.noexception.*;

class TemplateBuilder {
	FingerprintTransparency transparency = FingerprintTransparency.none;
	Cell size;
	Minutia[] minutiae;
	NeighborEdge[][] edges;
	void extract(byte[] image, double dpi) {
		// 지문 이미지에서 minutia를 추출하고 edge list구축

		// 모든 픽셀을 grayscale로 만든 1차원 double array를 만듬
		DoubleMap raw = readImage(image);
		if (Math.abs(dpi - 500) > Parameters.dpiTolerance)
			raw = scaleImage(raw, dpi);
		//transparency.logScaledImage(raw);
		// 이미지 크기
		size = raw.size();
		// 이미지크기를 개별 블록으로 쪼개어서 작업
		BlockMap blocks = new BlockMap(raw.width, raw.height, Parameters.blockSize);
		//transparency.logBlockMap(blocks);
		// 각 서브블록에 대해서 256개의 histogram을 구축한다.
		Histogram histogram = histogram(blocks, raw);
		// secondary size로 만들면서 인접 4개의 블록 histogram을 합한다.
		Histogram smoothHistogram = smoothHistogram(blocks, histogram);

		// 지문이 존재하는 영역을 표시
		BooleanMap mask = mask(blocks, histogram);

		// 해당 서브블록에 대해 주위블록들과 관계를 스무스하게 하는 것같음
		// 지문이 없는블록은 -1로 고정??
		DoubleMap equalized = equalize(blocks, raw, smoothHistogram, mask);

		// 모든 서브블록에 대해 진해지는 방향의 각도를 계산한다.
		DoubleMap orientation = orientationMap(equalized, mask, blocks);

		// 리지의 경계선에 수직방향으로 평균을낸다 (해상도32)
		// 진행방향으로 뭉개진다
		Cell[][] smoothedLines = orientedLines(Parameters.parallelSmoothinigResolution, Parameters.parallelSmoothinigRadius, Parameters.parallelSmoothinigStep);
		DoubleMap smoothed = smoothRidges(equalized, orientation, mask, blocks, 0, smoothedLines);
		//transparency.logParallelSmoothing(smoothed);

		// 리지의 진행방향으로 평균을 낸다 (해상도11)
		// 수직방향으로 뭉개진다
		Cell[][] orthogonalLines = orientedLines(Parameters.orthogonalSmoothinigResolution, Parameters.orthogonalSmoothinigRadius, Parameters.orthogonalSmoothinigStep);
		DoubleMap orthogonal = smoothRidges(smoothed, orientation, mask, blocks, Math.PI, orthogonalLines);

		//transparency.logOrthogonalSmoothing(orthogonal);
		// 바이너리 이미지로 만듬
		BooleanMap binary = binarize(smoothed, orthogonal, mask, blocks);
		// 지문이 있는 영역을 true로 채움
		BooleanMap pixelMask = fillBlocks(mask, blocks);
		// 바이너리 이미지의 지저분한 부분을 정리한다.
		cleanupBinarized(binary, pixelMask);
		//transparency.logPixelMask(pixelMask);
		// 지문 부분의 이미지를 역전한다.
		BooleanMap inverted = invert(binary, pixelMask);
		// inner mask를 기존 매스크를 줄여서 만든다??
		BooleanMap innerMask = innerMask(pixelMask);


		Skeleton ridges = new Skeleton(binary, SkeletonType.RIDGES, transparency);
		Skeleton valleys = new Skeleton(inverted, SkeletonType.VALLEYS, transparency);

		collectMinutiae(ridges, MinutiaType.ENDING);
		collectMinutiae(valleys, MinutiaType.BIFURCATION);

		//transparency.logSkeletonMinutiae(this);
		maskMinutiae(innerMask);
		removeMinutiaClouds();
		limitTemplateSize();
		shuffleMinutiae();
		buildEdgeTable();
	}
	void deserialize(String json) {
		JsonTemplate data = new Gson().fromJson(json, JsonTemplate.class);
		size = data.size();
		minutiae = data.minutiae();
		transparency.logDeserializedMinutiae(this);
		buildEdgeTable();
	}
	void convert(byte[] iso) {
		if (iso.length < 30)
			throw new IllegalArgumentException("Array too small to be an ISO 19794-2 template");
		try {
			DataInput in = new DataInputStream(new ByteArrayInputStream(iso));
			// 4B magic header "FMR\0"
			if (in.readByte() != 'F' || in.readByte() != 'M' || in.readByte() != 'R' || in.readByte() != 0)
				throw new IllegalArgumentException("This is not an ISO 19794-2 template");
			// 4B version " 20\0"
			// 4B template length in bytes (should be 28 + 6 * count + 2 + extra-data)
			// 2B junk
			in.skipBytes(10);
			// image size
			int width = in.readUnsignedShort();
			int height = in.readUnsignedShort();
			// pixels per cm X and Y, assuming 500dpi
			int xPixelsPerCM = in.readShort();
			int yPixelsPerCM = in.readShort();
			transparency.logIsoMetadata(width, height, xPixelsPerCM, yPixelsPerCM);
			double dpiX = xPixelsPerCM * 2.55;
			double dpiY = yPixelsPerCM * 2.55;
			boolean rescaleX = Math.abs(dpiX - 500) > Parameters.dpiTolerance;
			boolean rescaleY = Math.abs(dpiY - 500) > Parameters.dpiTolerance;
			if (rescaleX)
				width = (int)Math.round(width / dpiX * 500);
			if (rescaleY)
				height = (int)Math.round(height / dpiY * 500);
			size = new Cell(width, height);
			// 1B number of fingerprints in the template (assuming 1)
			// 1B junk
			// 1B finger position
			// 1B junk
			// 1B fingerprint quality
			in.skipBytes(5);
			// minutia count
			int count = in.readUnsignedByte();
			List<Minutia> list = new ArrayList<>();
			for (int i = 0; i < count; ++i) {
				// X position, upper two bits are type
				int packedX = in.readUnsignedShort();
				// Y position, upper two bits ignored
				int packedY = in.readUnsignedShort();
				// angle, 0..255 equivalent to 0..2pi
				int angle = in.readUnsignedByte();
				// 1B minutia quality
				in.skipBytes(1);
				// type: 01 ending, 10 bifurcation, 00 other (treated as ending)
				int type = (packedX >> 14) & 0x3;
				int x = packedX & 0x3fff;
				int y = packedY & 0x3fff;
				if (rescaleX)
					x = (int)Math.round(x / dpiX * 500);
				if (rescaleY)
					y = (int)Math.round(y / dpiY * 500);
				Minutia minutia = new Minutia(
					new Cell(x, y),
					Angle.complementary(angle * Angle.PI2 / 256.0),
					type == 2 ? MinutiaType.BIFURCATION : MinutiaType.ENDING);
				list.add(minutia);
			}
			// extra data length
			int extra = in.readUnsignedShort();
			// variable-length extra data section
			in.skipBytes(extra);
			minutiae = list.stream().toArray(Minutia[]::new);
			transparency.logIsoMinutiae(this);
		} catch (IOException e) {
			throw new IllegalArgumentException("Invalid ISO 19794-2 template", e);
		}
		shuffleMinutiae();
		buildEdgeTable();
	}
	DoubleMap readImage(byte[] serialized) {
		BufferedImage buffered = Exceptions.sneak().get(() -> ImageIO.read(new ByteArrayInputStream(serialized)));
		if (buffered == null)
			throw new IllegalArgumentException("Unsupported image format");
		// 이미지의 크기를 가져온다 
		int width = buffered.getWidth();
		int height = buffered.getHeight();
		// 이미지 크기만큼 저장공간확보
		int[] pixels = new int[width * height];
		// pixels로 데이터를 읽어냄
		buffered.getRGB(0, 0, width, height, pixels, 0, width);
		// double형식의 공간 준비
		DoubleMap map = new DoubleMap(width, height);
		// 모든 픽셀에 대해서 RGB를 분리한 후 grayscale image를 생성(그냥 평균???)
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				int pixel = pixels[y * width + x];
				int color = (pixel & 0xff) + ((pixel >> 8) & 0xff) + ((pixel >> 16) & 0xff);
				map.set(x, y, 1 - color * (1.0 / (3.0 * 255.0)));
			}
		}
		transparency.logDecodedImage(map);
		return map;
	}
	private DoubleMap scaleImage(DoubleMap input, double dpi) {
		return scaleImage(input, (int)Math.round(500.0 / dpi * input.width), (int)Math.round(500.0 / dpi * input.height));
	}
	private DoubleMap scaleImage(DoubleMap input, int newWidth, int newHeight) {
		// 이미지를 resize함
		DoubleMap output = new DoubleMap(newWidth, newHeight);
		double scaleX = newWidth / (double)input.width;
		double scaleY = newHeight / (double)input.height;
		double descaleX = 1 / scaleX;
		double descaleY = 1 / scaleY;
		for (int y = 0; y < newHeight; ++y) {
			double y1 = y * descaleY;
			double y2 = y1 + descaleY;
			int y1i = (int)y1;
			int y2i = (int)Math.ceil(y2);
			for (int x = 0; x < newWidth; ++x) {
				double x1 = x * descaleX;
				double x2 = x1 + descaleX;
				int x1i = (int)x1;
				int x2i = (int)Math.ceil(x2);
				double sum = 0;
				// 새 픽셀에 둘러싸인 모든 기존 픽셀을
				// 새 픽셀과의 거리를 고려하여 누적함
				for (int oy = y1i; oy < y2i; ++oy) {
					double ry = Math.min(oy + 1, y2) - Math.max(oy, y1);
					for (int ox = x1i; ox < x2i; ++ox) {
						double rx = Math.min(ox + 1, x2) - Math.max(ox, x1);
						sum += rx * ry * input.get(ox, oy);
					}
				}
				// 새 픽셀을 저장
				output.set(x, y, sum * (scaleX * scaleY));
			}
		}
		return output;
	}
	private Histogram histogram(BlockMap blocks, DoubleMap image) {
		// 모든 서브블록 별로 256개의 histogram을 구축한다.
		Histogram histogram = new Histogram(blocks.primary.blocks, Parameters.histogramDepth);
		for (Cell block : blocks.primary.blocks) {
			// 모든 15x15정도의 블록에 대해서 area block을 얻고
			Block area = blocks.primary.block(block);
			for (int y = area.top(); y < area.bottom(); ++y)
				for (int x = area.left(); x < area.right(); ++x) {
					// image가 0.0~1.0이기 때문에 depth 256를 곱하여
					// 정수화하면 0~255의 값이 나옴
					int depth = (int)(image.get(x, y) * histogram.depth);
					// 해당 블록의 histogram ++
					histogram.increment(block, histogram.constrain(depth));
				}
		}
		transparency.logHistogram(histogram);
		return histogram;
	}
	private Histogram smoothHistogram(BlockMap blocks, Histogram input) {
		Cell[] blocksAround = new Cell[] { new Cell(0, 0), new Cell(-1, 0), new Cell(0, -1), new Cell(-1, -1) };
		Histogram output = new Histogram(blocks.secondary.blocks, input.depth);
		for (Cell corner : blocks.secondary.blocks) {
			for (Cell relative : blocksAround) {
				Cell block = corner.plus(relative);
				if (blocks.primary.blocks.contains(block)) {
					for (int i = 0; i < input.depth; ++i)
						output.add(corner, i, input.get(block, i));
				}
			}
		}
		transparency.logSmoothedHistogram(output);
		return output;
	}
	private BooleanMap mask(BlockMap blocks, Histogram histogram) {
		// 블록히스토그램에서 상하위 8%씩을 제외한 영역의 비율계산
		// 넓을 수록 선명한 그림이란 의미임
		DoubleMap contrast = clipContrast(blocks, histogram);
		// 사용가능한 contrast영역이
		// 전체에서 17/255보다 작으면 그 블록을 체크한다.
		BooleanMap mask = filterAbsoluteContrast(contrast);
		// contrast가 평균에서 많이 떨어지는 것들도 체크한다.
		mask.merge(filterRelativeContrast(contrast, blocks));
		//transparency.logCombinedMask(mask);
		// r(9)영역안에 불량블록수가 일정 수 이상이면 체크
		mask.merge(vote(mask, null, Parameters.contrastVoteRadius, Parameters.contrastVoteMajority, Parameters.contrastVoteBorderDistance));
		// 해당 서브를록을 기준으로 4x4개의 블록중 70%가 불량이면 체크
		// 지문 외곽선을 따라 일부가 지문영역로 포함됨
		mask.merge(filterBlockErrors(mask));
		// 지문존재로 변경
		mask.invert();
		// 일부 외곽선을 따라 지문영역이 축소됨
		// 해당 서브를록을 기준으로 4x4개의 블록중 70%가 지문이 존재하면 체크
		mask.merge(filterBlockErrors(mask));
		// 해당 서브를록을 기준으로 4x4개의 블록중 70%가 지문이 존재하면 체크
		mask.merge(filterBlockErrors(mask));
		// r(9)영역안에 86%양품이면 양품
		mask.merge(vote(mask, null, Parameters.maskVoteRadius, Parameters.maskVoteMajority, Parameters.maskVoteBorderDistance));
		//transparency.logFilteredMask(mask);
		return mask;
	}
	private DoubleMap clipContrast(BlockMap blocks, Histogram histogram) {
		// 그리드 개수에 해당하는 result공간할당
		// 모든 블록에서 하위,상위 8%정도를 제외한 histogram영역의 비율을 계산
		DoubleMap result = new DoubleMap(blocks.primary.blocks);
		for (Cell block : blocks.primary.blocks) {
			// 모든 서브블록에 대해서
			// 공간 내 모든 histogram을 합하고
			int volume = histogram.sum(block);
			// 어두운 쪽에서 클립 level을 정하고 그 이하는 모두 black으로 바꿈
			// 밝은 쪽은 클립리미트 이상은 모두 white로 바꿈
			// 현재는 아래쪽 8%를 클립리미트로 설정
			int clipLimit = (int)Math.round(volume * Parameters.clippedContrast);
			int accumulator = 0;
			int lowerBound = histogram.depth - 1;
			for (int i = 0; i < histogram.depth; ++i) {
				accumulator += histogram.get(block, i);
				if (accumulator > clipLimit) {
					lowerBound = i;
					break;
				}
			}
			accumulator = 0;
			int upperBound = 0;
			for (int i = histogram.depth - 1; i >= 0; --i) {
				accumulator += histogram.get(block, i);
				if (accumulator > clipLimit) {
					upperBound = i;
					break;
				}
			}
			result.set(block, (upperBound - lowerBound) * (1.0 / (histogram.depth - 1)));
		}
		transparency.logClippedContrast(result);
		return result;
	}
	private BooleanMap filterAbsoluteContrast(DoubleMap contrast) {
		// 사용가능한 contrast영역이
		// 전체에서 17/255보다 작으면 그 블록을 체크한다.
		BooleanMap result = new BooleanMap(contrast.size());
		// 모든 contrast값에 대해
		for (Cell block : contrast.size())
			if (contrast.get(block) < Parameters.minAbsoluteContrast)
				result.set(block, true);
		transparency.logAbsoluteContrastMask(result);
		return result;
	}
	private BooleanMap filterRelativeContrast(DoubleMap contrast, BlockMap blocks) {
		// 서브블록 개수만큼의 배열에 모든 contrast값을 넣은 후
		List<Double> sortedContrast = new ArrayList<>();
		for (Cell block : contrast.size())
			sortedContrast.add(contrast.get(block));
		// 큰값이 앞으로 오도록 정렬
		sortedContrast.sort(Comparator.<Double>naturalOrder().reversed());
		// 서브블록당 평균적인 픽셀개수
		int pixelsPerBlock = blocks.pixels.area() / blocks.primary.blocks.area();
		// 410x410보다 큰 이미지이면 410x410정도의 개수만 취급
		int sampleCount = Math.min(sortedContrast.size(), Parameters.relativeContrastSample / pixelsPerBlock);
		// 
		int consideredBlocks = Math.max((int)Math.round(sampleCount * Parameters.relativeContrastPercentile), 1);
		// 모든 contrast의 평균을 구함 (limiting하여서)
		double averageContrast = sortedContrast.stream().mapToDouble(n -> n).limit(consideredBlocks).average().getAsDouble();
		// 평균의 아래쪽 34%는 limit로 지정
		double limit = averageContrast * Parameters.minRelativeContrast;
		// contrast가 limit이하인 것의 맵을 구축
		BooleanMap result = new BooleanMap(blocks.primary.blocks);
		for (Cell block : blocks.primary.blocks)
			if (contrast.get(block) < limit)
				result.set(block, true);
		transparency.logRelativeContrastMask(result);
		return result;
	}
	private BooleanMap vote(BooleanMap input, BooleanMap mask, int radius, double majority, int borderDistance) {
		//영역안에 불량서브블록 수가 일정 수 이상이면 체크한다.

		Cell size = input.size();
		// input에서 바깥쪽 7개씩을 뺀 영역을 잡음
		Block rect = new Block(borderDistance, borderDistance, size.x - 2 * borderDistance, size.y - 2 * borderDistance);
		// radius제곱만큼에서 여유를 두고 배열을 만듬
		int[] thresholds = IntStream.range(0, Integers.sq(2 * radius + 1) + 1).map(i -> (int)Math.ceil(majority * i)).toArray();
		// input(서브블록개수) 만큼의 intmap
		IntMap counts = new IntMap(size);
		// input(서브블록개수) 만큼
		BooleanMap output = new BooleanMap(size);
		for (int y = rect.top(); y < rect.bottom(); ++y) {
			// y를 기준으로 반경 radius를 영역을 잡고
			// 그림을 벗어나는 부분은 제외한다.
			int superTop = y - radius - 1;
			int superBottom = y + radius;
			int yMin = Math.max(0, y - radius);
			int yMax = Math.min(size.y - 1, y + radius);
			int yRange = yMax - yMin + 1;
			for (int x = rect.left(); x < rect.right(); ++x)
				// mask가 없거나 
				// 해당서브블록이 불량하면 
				if (mask == null || mask.get(x, y)) {
					// 좌측, 위쪽, 대각선 위쪽 값
					int left = x > 0 ? counts.get(x - 1, y) : 0;
					int top = y > 0 ? counts.get(x, y - 1) : 0;
					int diagonal = x > 0 && y > 0 ? counts.get(x - 1, y - 1) : 0;
					// x를 기준으로 반경영역을 잡고
					// 그림을 벗어나는 부분은 제외한다.
					int xMin = Math.max(0, x - radius);
					int xMax = Math.min(size.x - 1, x + radius);

					int ones;
					if (left > 0 && top > 0 && diagonal > 0) {
						ones = top + left - diagonal - 1;
						int superLeft = x - radius - 1;
						int superRight = x + radius;
						if (superLeft >= 0 && superTop >= 0 && input.get(superLeft, superTop))
							++ones;
						if (superLeft >= 0 && superBottom < size.y && input.get(superLeft, superBottom))
							--ones;
						if (superRight < size.x && superTop >= 0 && input.get(superRight, superTop))
							--ones;
						if (superRight < size.x && superBottom < size.y && input.get(superRight, superBottom))
							++ones;
					} else {
						// 왼쪽,위쪽 처음 라인이면
						// 반경 안의 모든 불량서브블록의 수를 센다
						ones = 0;
						for (int ny = yMin; ny <= yMax; ++ny)
							for (int nx = xMin; nx <= xMax; ++nx)
								if (input.get(nx, ny))
									++ones;
					}
					counts.set(x, y, ones + 1);
					// 불량 서브블록 수가 일정 수보다 크면 체크
					if (ones >= thresholds[yRange * (xMax - xMin + 1)])
						output.set(x, y, true);
				}
		}
		return output;
	}
	private BooleanMap filterBlockErrors(BooleanMap input) {
		return vote(input, null, Parameters.blockErrorsVoteRadius, Parameters.blockErrorsVoteMajority, Parameters.blockErrorsVoteBorderDistance);
	}
	private DoubleMap equalize(BlockMap blocks, DoubleMap image, Histogram histogram, BooleanMap blockMask) {
		final double rangeMin = -1;
		final double rangeMax = 1;
		final double rangeSize = rangeMax - rangeMin;
		final double widthMax = rangeSize / 256 * Parameters.maxEqualizationScaling;
		final double widthMin = rangeSize / 256 * Parameters.minEqualizationScaling;
		double[] limitedMin = new double[histogram.depth];
		double[] limitedMax = new double[histogram.depth];
		double[] dequantized = new double[histogram.depth];
		for (int i = 0; i < histogram.depth; ++i) {
			limitedMin[i] = Math.max(i * widthMin + rangeMin, rangeMax - (histogram.depth - 1 - i) * widthMax);
			limitedMax[i] = Math.min(i * widthMax + rangeMin, rangeMax - (histogram.depth - 1 - i) * widthMin);
			dequantized[i] = i / (double)(histogram.depth - 1);
		}
		Map<Cell, double[]> mappings = new HashMap<>();
		for (Cell corner : blocks.secondary.blocks) {
			double[] mapping = new double[histogram.depth];
			mappings.put(corner, mapping);
			// 현재 블록이나 이전, 위, 대각선위의 블록중에 하나라도 지문이있으면
			if (blockMask.get(corner, false) || blockMask.get(corner.x - 1, corner.y, false)
				|| blockMask.get(corner.x, corner.y - 1, false) || blockMask.get(corner.x - 1, corner.y - 1, false)) {
				double step = rangeSize / histogram.sum(corner);
				double top = rangeMin;
				for (int i = 0; i < histogram.depth; ++i) {
					double band = histogram.get(corner, i) * step;
					double equalized = top + dequantized[i] * band;
					top += band;
					if (equalized < limitedMin[i])
						equalized = limitedMin[i];
					if (equalized > limitedMax[i])
						equalized = limitedMax[i];
					mapping[i] = equalized;
				}
			}
		}
		DoubleMap result = new DoubleMap(blocks.pixels);
		for (Cell block : blocks.primary.blocks) {
			Block area = blocks.primary.block(block);
			if (blockMask.get(block)) {
				// 해당블록에 지문이 있으면
				// 오른쪽, 아래쪽, 대각선아래쪽의 테이블을 얻는다
				double[] topleft = mappings.get(block);
				double[] topright = mappings.get(new Cell(block.x + 1, block.y));
				double[] bottomleft = mappings.get(new Cell(block.x, block.y + 1));
				double[] bottomright = mappings.get(new Cell(block.x + 1, block.y + 1));
				for (int y = area.top(); y < area.bottom(); ++y)
					for (int x = area.left(); x < area.right(); ++x) {
						int depth = histogram.constrain((int)(image.get(x, y) * histogram.depth));
						double rx = (x - area.x + 0.5) / area.width;
						double ry = (y - area.y + 0.5) / area.height;
						result.set(x, y, Doubles.interpolate(bottomleft[depth], bottomright[depth], topleft[depth], topright[depth], rx, ry));
					}
			} else {
				for (int y = area.top(); y < area.bottom(); ++y)
					for (int x = area.left(); x < area.right(); ++x)
						result.set(x, y, -1);
			}
		}
		transparency.logEqualizedImage(result);
		return result;
	}
	private DoubleMap orientationMap(DoubleMap image, BooleanMap mask, BlockMap blocks) {
		// 모든 픽셀에 대해 해당 픽셀에서 진해지는 쪽 방향 벡터를 생성
		PointMap accumulated = pixelwiseOrientation(image, mask, blocks);
		// 서브블록 별로 방향벡터를 누적한다.
		PointMap byBlock = blockOrientations(accumulated, blocks, mask);
		// 하나 이웃인 모든 블록의 방향벡터를 누적한다.
		PointMap smooth = smoothOrientation(byBlock, mask);
		// 모든 서브블록의 방향벡터를 각도로 환산한다.
		// 지문이 없는 영역은 0
		return orientationAngles(smooth, mask);
	}
	private static class ConsideredOrientation {
		// 정수 벡터
		Cell offset;
		// 해당 정수벡터 방향의 유닛벡터
		Point orientation;
	}
	private static class OrientationRandom {
		// 30bit 유사 난수 생성
		static final int prime = 1610612741;
		static final int bits = 30;
		static final int mask = (1 << bits) - 1;
		static final double scaling = 1.0 / (1 << bits);
		long state = prime * prime * prime;
		double next() {
			state *= prime;
			return ((state & mask) + 0.5) * scaling;
		}
	}
	private ConsideredOrientation[][] planOrientations() {
		OrientationRandom random = new OrientationRandom();
		// 50x20개의 난수(각도,거리)
		ConsideredOrientation[][] splits = new ConsideredOrientation[Parameters.orientationSplit][];
		for (int i = 0; i < Parameters.orientationSplit; ++i) {
			ConsideredOrientation[] orientations = splits[i] = new ConsideredOrientation[Parameters.orientationsChecked];
			for (int j = 0; j < Parameters.orientationsChecked; ++j) {
				ConsideredOrientation sample = orientations[j] = new ConsideredOrientation();
				do {
					// 0~180도 사이의 난수 각도 
					double angle = random.next() * Math.PI;
					// 2~6사이의 난수 거리
					double distance = Doubles.interpolateExponential(Parameters.minOrientationRadius, Parameters.maxOrientationRadius, random.next());
					// 해당 각도 거리의 위치
					sample.offset = Angle.toVector(angle).multiply(distance).round();
					// 0이나 중복이 없도록 반복 
				} while (sample.offset.equals(Cell.zero) || sample.offset.y < 0 || Arrays.stream(orientations).limit(j).anyMatch(o -> o.offset.equals(sample.offset)));
				sample.orientation = Angle.toVector(
						Angle.add(
							Angle.toOrientation(
								Angle.atan(sample.offset.toPoint())	// -pi/2 ~ pi/2
							),										// -pi ~ pi
						   	Math.PI
						)											// 0 ~ 2*pi
					);
			}
		}
		return splits;
	}
	private PointMap pixelwiseOrientation(DoubleMap input, BooleanMap mask, BlockMap blocks) {
		// 50x20개의 난수방향 생성
		ConsideredOrientation[][] neighbors = planOrientations();
		// 모든 픽셀에 대해 방향을 가진 맵 생성
		PointMap orientation = new PointMap(input.size());
		// y방향 서브블록들에 대해
		for (int blockY = 0; blockY < blocks.primary.blocks.y; ++blockY) {
			// 해당 row에서 실제 지문이 존재하는 영역을 계산
			Range maskRange = maskRange(mask, blockY);
			// 영역이 존재하면(빈row는 skip)
			if (maskRange.length() > 0) {
				// 해당 row에서 지문이 존재하는 영역의 pixel단위 영역
				Range validXRange = new Range(
					blocks.primary.block(maskRange.start, blockY).left(),
					blocks.primary.block(maskRange.end - 1, blockY).right());
				// 해당 블록의 모든 y에 대해
				for (int y = blocks.primary.block(0, blockY).top(); y < blocks.primary.block(0, blockY).bottom(); ++y) {
					// 난수 방향 20개에 대해
					for (ConsideredOrientation neighbor : neighbors[y % neighbors.length]) {
						// 난수방향의 x,y중 큰것을 반지름으로 하여
						int radius = Math.max(Math.abs(neighbor.offset.x), Math.abs(neighbor.offset.y));
						// y +- radius가 이력 그림안에 포함되면
						if (y - radius >= 0 && y + radius < input.height) {
							// 반지름을 고려하여 valid range를 수축한다.
							// +radius ~ max-radius
							Range xRange = new Range(Math.max(radius, validXRange.start), Math.min(input.width - radius, validXRange.end));
							// 영역안의 모든 픽셀에 대해
							for (int x = xRange.start; x < xRange.end; ++x) {
								// 해당 방향의 픽셀과 반대방향의 픽셀중 큰놈이 
								// 현재 픽셀값보다 크면
								// 해당 방향에 차이값만큼 스케일해서 더한다.
								double before = input.get(x - neighbor.offset.x, y - neighbor.offset.y);
								double at = input.get(x, y);
								double after = input.get(x + neighbor.offset.x, y + neighbor.offset.y);
								double strength = at - Math.max(before, after);
								if (strength > 0)
									orientation.add(x, y, neighbor.orientation.multiply(strength));
							}
						}
					}
				}
			}
		}
		transparency.logPixelwiseOrientation(orientation);
		return orientation;
	}
	private static Range maskRange(BooleanMap mask, int y) {
		int first = -1;
		int last = -1;
		for (int x = 0; x < mask.width; ++x)
			if (mask.get(x, y)) {
				last = x;
				if (first < 0)
					first = x;
			}
		if (first >= 0)
			return new Range(first, last + 1);
		else
			return Range.zero;
	}
	private PointMap blockOrientations(PointMap orientation, BlockMap blocks, BooleanMap mask) {
		// 서브블록 크기의 맵을 생성
		PointMap sums = new PointMap(blocks.primary.blocks);
		// 모든 서브블록에 대해
		for (Cell block : blocks.primary.blocks) {
			// 지문이 존재하는 영역이면
			if (mask.get(block)) {
				// 해당 서브블록 영역을 얻고
				Block area = blocks.primary.block(block);
				// 그영역의 모든 픽셀의 방향을 벡터 합한다.
				for (int y = area.top(); y < area.bottom(); ++y)
					for (int x = area.left(); x < area.right(); ++x)
						sums.add(block, orientation.get(x, y));
			}
		}
		transparency.logBlockOrientation(sums);
		return sums;
	}
	private PointMap smoothOrientation(PointMap orientation, BooleanMap mask) {
		// size는 서브블록개수
		Cell size = mask.size();
		// 서브블록 개수만큼의
		PointMap smoothed = new PointMap(size);
		for (Cell block : size)
			// 지문이 존재하는 영역이면
			if (mask.get(block)) {
				// 해당 서브블록기준으로 +-1영역의 블록중 이미지 안에 존재하는 블록들
				Block neighbors = Block.around(block, Parameters.orientationSmoothingRadius).intersect(new Block(size));
				// 이웃블록의 모든 방향의 벡터합
				for (int ny = neighbors.top(); ny < neighbors.bottom(); ++ny)
					for (int nx = neighbors.left(); nx < neighbors.right(); ++nx)
						if (mask.get(nx, ny))
							smoothed.add(block, orientation.get(nx, ny));
			}
		transparency.logSmoothedOrientation(smoothed);
		return smoothed;
	}
	private static DoubleMap orientationAngles(PointMap vectors, BooleanMap mask) {
		Cell size = mask.size();
		DoubleMap angles = new DoubleMap(size);
		for (Cell block : size)
			if (mask.get(block))
				angles.set(block, Angle.atan(vectors.get(block)));
		return angles;
	}
	private Cell[][] orientedLines(int resolution, int radius, double step) {
		// resolution만큼 방향에 대해 radius부터 작아지는 방향으로 step씩 당겨가며 line을 만듬
		// cell들의 집합인데 점선이 되지 않을까??
		Cell[][] result = new Cell[resolution][];
		for (int orientationIndex = 0; orientationIndex < resolution; ++orientationIndex) {
			List<Cell> line = new ArrayList<>();
			line.add(Cell.zero);
			// 단위벡터
			Point direction = Angle.toVector(
					// 0~pi까지로 축소
					Angle.fromOrientation(
						// 0~ 2pi까지를 resoultion으로 등분한 각도
						Angle.bucketCenter(orientationIndex, resolution)
					)
				);
			// 반지름을 /step 하여 0.5가 될 때까지
			for (double r = radius; r >= 0.5; r /= step) {
				// 해당 방향으로 반지름반큼 진행한 위치
				Cell sample = direction.multiply(r).round();
				if (!line.contains(sample)) {
					line.add(sample);
					line.add(sample.negate());
				}
			}
			result[orientationIndex] = line.toArray(new Cell[line.size()]);
		}
		return result;
	}
	private static DoubleMap smoothRidges(DoubleMap input, DoubleMap orientation, BooleanMap mask, BlockMap blocks, double angle, Cell[][] lines) {
		// 모든 픽셀에 대해 output 공간을 할당
		DoubleMap output = new DoubleMap(input.size());
		// 모든 서브블록에 대해
		for (Cell block : blocks.primary.blocks) {
			// 지문이 존재하는 영역이면
			if (mask.get(block)) {
				// 해당 블록의 방향 + angle에 가장 근접한 라인 배열을 선택
				Cell[] line = lines[Angle.quantize(
										Angle.add(
											orientation.get(block), angle	// 해당 블록의 orientation에 각도를 더하고
										), lines.length				// 라인 해상도로 분리하여
									)];
				// 모든 라인 포인트에 대해
				for (Cell linePoint : line) {
					Block target = blocks.primary.block(block);
					// 라인의 점들을 해당 서브블록으로 옮기고 블록을 벗어난 부분 제거
					Block source = target.move(linePoint).intersect(new Block(blocks.pixels));
					// 원래 위치로 되돌리고
					target = source.move(linePoint.negate());
					// 리지의 경계선의 수직방향의 여러 픽셀값을 더한다.
					for (int y = target.top(); y < target.bottom(); ++y)
						for (int x = target.left(); x < target.right(); ++x)
							output.add(x, y, input.get(x + linePoint.x, y + linePoint.y));
				}
				// 여러값을 더했으므로 다시 나누어서 0.0 ~ 1.0영역으로 바꾼다.
				Block blockArea = blocks.primary.block(block);
				for (int y = blockArea.top(); y < blockArea.bottom(); ++y)
					for (int x = blockArea.left(); x < blockArea.right(); ++x)
						output.multiply(x, y, 1.0 / line.length);
			}
		}
		return output;
	}
	private BooleanMap binarize(DoubleMap input, DoubleMap baseline, BooleanMap mask, BlockMap blocks) {
		Cell size = input.size();
		// 전체 이미지 크기의 공간 확보
		BooleanMap binarized = new BooleanMap(size);
		for (Cell block : blocks.primary.blocks)
			if (mask.get(block)) {
				Block rect = blocks.primary.block(block);
				for (int y = rect.top(); y < rect.bottom(); ++y)
					for (int x = rect.left(); x < rect.right(); ++x)
						if (input.get(x, y) - baseline.get(x, y) > 0)
							binarized.set(x, y, true);
			}
		transparency.logBinarizedImage(binarized);
		return binarized;
	}
	private void cleanupBinarized(BooleanMap binary, BooleanMap mask) {
		Cell size = binary.size();
		BooleanMap inverted = new BooleanMap(binary);
		inverted.invert();
		// 조그만 섬
		BooleanMap islands = vote(inverted, mask, Parameters.binarizedVoteRadius, Parameters.binarizedVoteMajority, Parameters.binarizedVoteBorderDistance);
		// 조그만 구멍
		BooleanMap holes = vote(binary, mask, Parameters.binarizedVoteRadius, Parameters.binarizedVoteMajority, Parameters.binarizedVoteBorderDistance);
		for (int y = 0; y < size.y; ++y)
			for (int x = 0; x < size.x; ++x)
				// 섬과 구멍을 제거
				binary.set(x, y, binary.get(x, y) && !islands.get(x, y) || holes.get(x, y));
		// 대각선으로 다른 색깔인 경우 모두다 false로 만든다.
		removeCrosses(binary);
		transparency.logFilteredBinarydImage(binary);
	}
	private static void removeCrosses(BooleanMap input) {
		Cell size = input.size();
		boolean any = true;
		while (any) {
			any = false;
			for (int y = 0; y < size.y - 1; ++y)
				for (int x = 0; x < size.x - 1; ++x)
					if (input.get(x, y) && input.get(x + 1, y + 1) && !input.get(x, y + 1) && !input.get(x + 1, y)
						|| input.get(x, y + 1) && input.get(x + 1, y) && !input.get(x, y) && !input.get(x + 1, y + 1)) {
						input.set(x, y, false);
						input.set(x, y + 1, false);
						input.set(x + 1, y, false);
						input.set(x + 1, y + 1, false);
						any = true;
					}
		}
	}
	private static BooleanMap fillBlocks(BooleanMap mask, BlockMap blocks) {
		BooleanMap pixelized = new BooleanMap(blocks.pixels);
		for (Cell block : blocks.primary.blocks)
			if (mask.get(block))
				for (Cell pixel : blocks.primary.block(block))
					pixelized.set(pixel, true);
		return pixelized;
	}
	private static BooleanMap invert(BooleanMap binary, BooleanMap mask) {
		Cell size = binary.size();
		BooleanMap inverted = new BooleanMap(size);
		for (int y = 0; y < size.y; ++y)
			for (int x = 0; x < size.x; ++x)
				inverted.set(x, y, !binary.get(x, y) && mask.get(x, y));
		return inverted;
	}
	private BooleanMap innerMask(BooleanMap outer) {
		Cell size = outer.size();
		BooleanMap inner = new BooleanMap(size);
		for (int y = 1; y < size.y - 1; ++y)
			for (int x = 1; x < size.x - 1; ++x)
				inner.set(x, y, outer.get(x, y));
		if (Parameters.innerMaskBorderDistance >= 1)
			inner = shrinkMask(inner, 1);
		int total = 1;
		for (int step = 1; total + step <= Parameters.innerMaskBorderDistance; step *= 2) {
			inner = shrinkMask(inner, step);
			total += step;
		}
		if (total < Parameters.innerMaskBorderDistance)
			inner = shrinkMask(inner, Parameters.innerMaskBorderDistance - total);
		transparency.logInnerMask(inner);
		return inner;
	}
	private static BooleanMap shrinkMask(BooleanMap mask, int amount) {
		Cell size = mask.size();
		BooleanMap shrunk = new BooleanMap(size);
		for (int y = amount; y < size.y - amount; ++y)
			for (int x = amount; x < size.x - amount; ++x)
				shrunk.set(x, y, mask.get(x, y - amount) && mask.get(x, y + amount) && mask.get(x - amount, y) && mask.get(x + amount, y));
		return shrunk;
	}
	private void collectMinutiae(Skeleton skeleton, MinutiaType type) {
		minutiae = Stream.concat(
			Arrays.stream(Optional.ofNullable(minutiae).orElse(new Minutia[0])),
			skeleton.minutiae.stream()
				.filter(m -> m.ridges.size() == 1)
				.map(m -> new Minutia(m.position, m.ridges.get(0).direction(), type)))
			.toArray(Minutia[]::new);
	}
	private void maskMinutiae(BooleanMap mask) {
		minutiae = Arrays.stream(minutiae)
			.filter(minutia -> {
				Cell arrow = Angle.toVector(minutia.direction).multiply(-Parameters.maskDisplacement).round();
				return mask.get(minutia.position.plus(arrow), false);
			})
			.toArray(Minutia[]::new);
		transparency.logInnerMinutiae(this);
	}
	private void removeMinutiaClouds() {
		int radiusSq = Integers.sq(Parameters.minutiaCloudRadius);
		Set<Minutia> removed = Arrays.stream(minutiae)
			.filter(minutia -> Parameters.maxCloudSize < Arrays.stream(minutiae)
				.filter(neighbor -> neighbor.position.minus(minutia.position).lengthSq() <= radiusSq)
				.count() - 1)
			.collect(toSet());
		minutiae = Arrays.stream(minutiae)
			.filter(minutia -> !removed.contains(minutia))
			.toArray(Minutia[]::new);
		transparency.logRemovedMinutiaClouds(this);
	}
	private void limitTemplateSize() {
		if (minutiae.length > Parameters.maxMinutiae) {
			minutiae = Arrays.stream(minutiae)
				.sorted(Comparator.<Minutia>comparingInt(
					minutia -> Arrays.stream(minutiae)
						.mapToInt(neighbor -> minutia.position.minus(neighbor.position).lengthSq())
						.sorted()
						.skip(Parameters.sortByNeighbor)
						.findFirst().orElse(Integer.MAX_VALUE))
					.reversed())
				.limit(Parameters.maxMinutiae)
				.toArray(Minutia[]::new);
		}
		transparency.logTopMinutiae(this);
	}
	private void shuffleMinutiae() {
		int prime = 1610612741;
		Arrays.sort(minutiae, Comparator
			.comparing((Minutia m) -> ((m.position.x * prime) + m.position.y) * prime)
			.thenComparing(m -> m.position.x)
			.thenComparing(m -> m.position.y)
			.thenComparing(m -> m.direction)
			.thenComparing(m -> m.type));
		transparency.logShuffledMinutiae(this);
	}
	private void buildEdgeTable() {
		edges = new NeighborEdge[minutiae.length][];
		List<NeighborEdge> star = new ArrayList<>();
		int[] allSqDistances = new int[minutiae.length];
		// 모든 minutia에 대하여
		for (int reference = 0; reference < edges.length; ++reference) {
			//현재 minutia의 위치
			Cell referencePosition = minutiae[reference].position;
			// 490*490
			int sqMaxDistance = Integers.sq(Parameters.edgeTableRange);
			// 총 minutia의 개수가 9이상이면
			if (minutiae.length - 1 > Parameters.edgeTableNeighbors) {
				// 모든 minutia에 대하여
				// 거리의 제곱을 각각 계산
				for (int neighbor = 0; neighbor < minutiae.length; ++neighbor)
					allSqDistances[neighbor] = referencePosition.minus(minutiae[neighbor].position).lengthSq();
				// 정렬
				Arrays.sort(allSqDistances);
				// 가까운 순서로 9개째를 max로 선택
				sqMaxDistance = allSqDistances[Parameters.edgeTableNeighbors];
			}
			// 모든 minutia에대해서(reference 제외)
			// distance가 가까운쪽으로 10개까지만 선택 (0~9)
			for (int neighbor = 0; neighbor < minutiae.length; ++neighbor) {
				if (neighbor != reference && referencePosition.minus(minutiae[neighbor].position).lengthSq() <= sqMaxDistance)
					star.add(new NeighborEdge(minutiae, reference, neighbor));
			}
			// edge의 길이로 소팅하고 같으면 neighbor번호로 소팅
			star.sort(Comparator.<NeighborEdge>comparingInt(e -> e.length).thenComparingInt(e -> e.neighbor));
			// edge 리스트의 크기를 제한??
			while (star.size() > Parameters.edgeTableNeighbors)
				star.remove(star.size() - 1);
			//최종 edges는 Array형으로 변환
			edges[reference] = star.toArray(new NeighborEdge[star.size()]);
			star.clear();
		}
		transparency.logEdgeTable(edges);
	}
}
