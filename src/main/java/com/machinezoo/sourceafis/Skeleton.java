// Part of SourceAFIS: https://sourceafis.machinezoo.com
package com.machinezoo.sourceafis;

import java.nio.*;
import java.util.*;

class Skeleton {
	private final FingerprintTransparency logger;
	final SkeletonType type;
	final Cell size;
	final List<SkeletonMinutia> minutiae = new ArrayList<>();
	Skeleton(BooleanMap binary, SkeletonType type, FingerprintTransparency logger) {
		this.type = type;
		this.logger = logger;
		logger.logBinarizedSkeleton(type, binary);
		size = binary.size();
		// thin이미지 생성
		BooleanMap thinned = thin(binary);
		// thined이미지에서 모든 픽셀에 대해 이웃이 1이면 ending
		// 이웃이 2보다 많으면 bifurcation
		List<Cell> minutiaPoints = findMinutiae(thinned);
		// 이웃이 모두 미누셔일 때 이웃들을 엮는다
		// 아마도 ending-ending, ending-bifurcation-ending의 경우밖에 없지 않을까??
		Map<Cell, List<Cell>> linking = linkNeighboringMinutiae(minutiaPoints);
		// 미누셔의 이웃을 정리하여 독립적인 미누셔만 남긴다.
		Map<Cell, SkeletonMinutia> minutiaMap = minutiaCenters(linking);
		// 각각의 미누셔에서 시작되고 끝나는 모든 리지의 픽셀들을 추적하여 데이터베이스를 만든다
		traceRidges(thinned, minutiaMap);
		// 리지의 시작점과 미누셔가 일치하지 않는 경우 선을 긋는다.
		fixLinkingGaps();
		logger.logTracedSkeleton(this);
		//
		filter();
	}
	private enum NeighborhoodType {
		Skeleton,
		Ending,
		Removable
	}
	private BooleanMap thin(BooleanMap input) {
		// 이웃의 배치에 따라 현재 픽셀을 제거할 수 있는지를 미리 계산하여 둔다
		NeighborhoodType[] neighborhoodTypes = neighborhoodTypes();
		// 중간 값
		BooleanMap partial = new BooleanMap(size);
		for (int y = 1; y < size.y - 1; ++y)
			for (int x = 1; x < size.x - 1; ++x)
				partial.set(x, y, input.get(x, y));
		// 최종 결과
		// 초기값 아마도 false
		BooleanMap thinned = new BooleanMap(size);
		boolean removedAnything = true;
		// 여러번 반복 (현재 26번)
		// 더이상 제거할 게 없으면 중지
		for (int i = 0; i < Parameters.thinningIterations && removedAnything; ++i) {
			removedAnything = false;
			// 2x2영역에 대해
			for (int evenY = 0; evenY < 2; ++evenY)
				for (int evenX = 0; evenX < 2; ++evenX)
					// 전체 이미지를 2스텝식 이동하며 처리
					for (int y = 1 + evenY; y < size.y - 1; y += 2)
						for (int x = 1 + evenX; x < size.x - 1; x += 2)
							// 중간값에는 있고 최종 결과에 억는 픽셀이
							// 왼쪽만 픽셀이 없으면
							if (partial.get(x, y) && !thinned.get(x, y) && !(partial.get(x, y - 1) && partial.get(x, y + 1) && partial.get(x - 1, y) && partial.get(x + 1, y))) {
								// 주위 8개의 값을 비트맵으로
								int neighbors = (partial.get(x + 1, y + 1) ? 128 : 0)
									| (partial.get(x, y + 1) ? 64 : 0)
									| (partial.get(x - 1, y + 1) ? 32 : 0)
									| (partial.get(x + 1, y) ? 16 : 0)
									| (partial.get(x - 1, y) ? 8 : 0)
									| (partial.get(x + 1, y - 1) ? 4 : 0)
									| (partial.get(x, y - 1) ? 2 : 0)
									| (partial.get(x - 1, y - 1) ? 1 : 0);
								// 주위 8개의 픽셀로 미리 정보를 만들어 놓았다
								// 제거할 수 있으면 제거
								// ending의 경우 진짜 ending인지 한번더 체크
								if (neighborhoodTypes[neighbors] == NeighborhoodType.Removable
									|| neighborhoodTypes[neighbors] == NeighborhoodType.Ending
										&& isFalseEnding(partial, new Cell(x, y))) {
									removedAnything = true;
									partial.set(x, y, false);
								} else
									// 제거할 수 없는 경우이면 thinned에 업데이트
									thinned.set(x, y, true);
							}
		}
		logger.logThinnedSkeleton(type, thinned);
		return thinned;
	}
	private static NeighborhoodType[] neighborhoodTypes() {
		NeighborhoodType[] types = new NeighborhoodType[256];
		for (int mask = 0; mask < 256; ++mask) {
			boolean TL = (mask & 1) != 0;
			boolean TC = (mask & 2) != 0;
			boolean TR = (mask & 4) != 0;
			boolean CL = (mask & 8) != 0;
			boolean CR = (mask & 16) != 0;
			boolean BL = (mask & 32) != 0;
			boolean BC = (mask & 64) != 0;
			boolean BR = (mask & 128) != 0;
			int count = Integer.bitCount(mask);
			// 대각쪽 만 있는 경우
			boolean diagonal = !TC && !CL && TL || !CL && !BC && BL || !BC && !CR && BR || !CR && !TC && TR;
			// 좌우로 연결된 경우
			boolean horizontal = !TC && !BC && (TR || CR || BR) && (TL || CL || BL);
			// 아래위로 연결된 경우
			boolean vertical = !CL && !CR && (TL || TC || TR) && (BL || BC || BR);
			boolean end = (count == 1);
			if (end)
				types[mask] = NeighborhoodType.Ending;
			else if (!diagonal && !horizontal && !vertical)
				types[mask] = NeighborhoodType.Removable;
			else
				types[mask] = NeighborhoodType.Skeleton;
		}
		return types;
	}
	private static boolean isFalseEnding(BooleanMap binary, Cell ending) {
		for (Cell relativeNeighbor : Cell.cornerNeighbors) {
			Cell neighbor = ending.plus(relativeNeighbor);
			if (binary.get(neighbor)) {
				int count = 0;
				for (Cell relative2 : Cell.cornerNeighbors)
					if (binary.get(neighbor.plus(relative2), false))
						++count;
				return count > 2;
			}
		}
		return false;
	}
	private List<Cell> findMinutiae(BooleanMap thinned) {
		List<Cell> result = new ArrayList<>();
		for (Cell at : size)
			if (thinned.get(at)) {
				int count = 0;
				for (Cell relative : Cell.cornerNeighbors)
					if (thinned.get(at.plus(relative), false))
						++count;
				// 1이면 ending, 1이상이면 bifurcation
				if (count == 1 || count > 2)
					result.add(at);
			}
		return result;
	}
	// 이웃픽셀이 둘다 미누셔일 때 이웃으로 지정한다. 
	private static Map<Cell, List<Cell>> linkNeighboringMinutiae(List<Cell> minutiae) {
		Map<Cell, List<Cell>> linking = new HashMap<>();
		// 모든 미뉴셔에 대해
		for (Cell minutiaPos : minutiae) {
			List<Cell> ownLinks = null;
			// 해당 픽셀의 이웃픽셀들에 대해
			for (Cell neighborRelative : Cell.cornerNeighbors) {
				Cell neighborPos = minutiaPos.plus(neighborRelative);
				// 이웃 셀이 이미 미뉴셔맵에 있으면
				if (linking.containsKey(neighborPos)) {
					// 찾아진 미누셔의 이웃 리스트를 가져온다
					List<Cell> neighborLinks = linking.get(neighborPos);
					// 가져온 리스트가 현재리스트와 다르면
					// 여러방향으로 이웃일 때 맨처음 한번만 매치된다
					// 아니면 다른 쪽으로 이웃이 있거나
					if (neighborLinks != ownLinks) {
						// 이미 다른 쪽으로 이웃이 있을 때
						if (ownLinks != null) {
							// 현재까지 이웃을 모두 합치고
							neighborLinks.addAll(ownLinks);
							// 해당 이웃들에게 모두 새 이웃을 등록한다.
							for (Cell mergedPos : ownLinks)
								linking.put(mergedPos, neighborLinks);
						}
						ownLinks = neighborLinks;
					}
				}
			}
			// 처음에는 null이므로 새 리스트 할당 
			if (ownLinks == null)
				ownLinks = new ArrayList<>();
			// 현재 미누셔를 추가하고 
			ownLinks.add(minutiaPos);
			// 맵에 넣는다
			linking.put(minutiaPos, ownLinks);
		}
		return linking;
	}
	private Map<Cell, SkeletonMinutia> minutiaCenters(Map<Cell, List<Cell>> linking) {
		Map<Cell, SkeletonMinutia> centers = new HashMap<>();
		// 모든 미뉴셔에 대해
		for (Cell currentPos : linking.keySet()) {
			// 이웃리스트를 가져온다
			List<Cell> linkedMinutiae = linking.get(currentPos);
			// 맨 앞의 것을 primary로 지정
			Cell primaryPos = linkedMinutiae.get(0);
			// primary가 center에 들어 있지 않으면
			if (!centers.containsKey(primaryPos)) {
				// 
				Cell sum = Cell.zero;
				// 현재 미누셔와 이웃인 모든 미누셔에 대해
				// 모든 위치를 벡터섬한다.
				for (Cell linkedPos : linkedMinutiae)
					sum = sum.plus(linkedPos);
				// 모든 위치의 평균위치를 센터로 지정한다.
				Cell center = new Cell(sum.x / linkedMinutiae.size(), sum.y / linkedMinutiae.size());
				// SkeleconMinutia를 만든다??
				SkeletonMinutia minutia = new SkeletonMinutia(center);
				// 미 미누셔를 확정한다
				addMinutia(minutia);
				// 센터 맵에 넣고 중복처리를 하지 않도록 한다.
				centers.put(primaryPos, minutia);
			}
			//현재 위치도 centers에 넣어 중복처리를 막는다
			centers.put(currentPos, centers.get(primaryPos));
		}
		return centers;
	}
	private void traceRidges(BooleanMap thinned, Map<Cell, SkeletonMinutia> minutiaePoints) {
		//
		Map<Cell, SkeletonRidge> leads = new HashMap<>();
		// 모든 미누셔 포인트에 대해
		for (Cell minutiaPoint : minutiaePoints.keySet()) {
			// 이웃 8개의 방향에 대해
			for (Cell startRelative : Cell.cornerNeighbors) {
				Cell start = minutiaPoint.plus(startRelative);
				// 해당 이웃이 리지이고
				// 미누셔가 아니고(이건 당연)
				// 이미 처리되지 않았으면
				if (thinned.get(start, false) && !minutiaePoints.containsKey(start) && !leads.containsKey(start)) {
					// 리지에 현재위치와 찾아진 이웅을 추가한다.
					SkeletonRidge ridge = new SkeletonRidge();
					ridge.points.add(minutiaPoint);
					ridge.points.add(start);
					// 리지를 따라가며 계속 리지를 추적한다.
					Cell previous = minutiaPoint;
					Cell current = start;
					do {
						Cell next = Cell.zero;
						for (Cell nextRelative : Cell.cornerNeighbors) {
							next = current.plus(nextRelative);
							// 새 리치 픽셀이 찾아지면
							if (thinned.get(next, false) && !next.equals(previous))
								break;
						}
						previous = current;
						current = next;
						// 리지 리스트에 추가하고
						ridge.points.add(current);
						// 다른 리지가 나올 때까지 리지를 계속 추적한다.
					} while (!minutiaePoints.containsKey(current));
					// 찾아진 반대편 미누셔
					Cell end = current;
					// 리지에 시작점과 끝점을 추가하고
					ridge.start(minutiaePoints.get(minutiaPoint));
					ridge.end(minutiaePoints.get(end));
					// lead에 리지의 각 끝점에 대해 리지를 추가한다.
					leads.put(ridge.points.get(1), ridge);
					leads.put(ridge.reversed.points.get(1), ridge);
				}
			}
		}
	}
	private void fixLinkingGaps() {
		// 모든 미누셔에 대해
		for (SkeletonMinutia minutia : minutiae) {
			// 미누셔에 연결된 모든 리지에 대해
			for (SkeletonRidge ridge : minutia.ridges) {
				// 만약 리지의 첫번째 포인트가 현재 미누셔가 아니면
				if (!ridge.points.get(0).equals(minutia.position)) {
					// 현재 포인트부터 리지 시작점까지 연결선을 그린다??
					Cell[] filling = ridge.points.get(0).lineTo(minutia.position);
					for (int i = 1; i < filling.length; ++i)
						ridge.reversed.points.add(filling[i]);
				}
			}
		}
	}
	private void filter() {
		// 리지가 없는 모든 미누셔 제거
		removeDots();
		logger.logRemovedDots(this);
		// 서로다른 리지가 같은 양끝을 가지고 짧을 때 기공으로 간주하고 제거한다.
		removePores();
		removeGaps();
		removeTails();
		removeFragments();
	}
	private void removeDots() {
		List<SkeletonMinutia> removed = new ArrayList<>();
		// 현재 모든 미누셔에 대해
		for (SkeletonMinutia minutia : minutiae)
			// 리지가 없는 미누셔는 removed에 등록
			if (minutia.ridges.isEmpty())
				removed.add(minutia);
		// removed 모든 미누셔 제거
		for (SkeletonMinutia minutia : removed)
			removeMinutia(minutia);
	}
	private void removePores() {
		// 기공을 제거함
		// 모든 미누셔에 대해
		for (SkeletonMinutia minutia : minutiae) {
			// 리지가 3개 연결된 경우
			if (minutia.ridges.size() == 3) {
				// 모든 리지에 대해
				for (int exit = 0; exit < 3; ++exit) {
					// 현재 리지
					SkeletonRidge exitRidge = minutia.ridges.get(exit);
					// 다음 리지 2개
					SkeletonRidge arm1 = minutia.ridges.get((exit + 1) % 3);
					SkeletonRidge arm2 = minutia.ridges.get((exit + 2) % 3);
					// 현재 리지 외에 두개의 리지가 같은 끝을 가질 때
					// 그리고 리지가 되돌아오는 경우가 아닐 때
					if (arm1.end() == arm2.end() && exitRidge.end() != arm1.end() && arm1.end() != minutia && exitRidge.end() != minutia) {
						SkeletonMinutia end = arm1.end();
						// 중복된 상대편도 가지가 3개이고 리지의 길이가 짧을 때
						// 기공으로 간주하고 제거하고 리지를 서로 연결한다.
						// 미누셔 제거는 어디서??
						if (end.ridges.size() == 3 && arm1.points.size() <= Parameters.maxPoreArm && arm2.points.size() <= Parameters.maxPoreArm) {
							arm1.detach();
							arm2.detach();
							SkeletonRidge merged = new SkeletonRidge();
							merged.start(minutia);
							merged.end(end);
							for (Cell point : minutia.position.lineTo(end.position))
								merged.points.add(point);
						}
						break;
					}
				}
			}
		}
		removeKnots();
		logger.logRemovedPores(this);
	}
	private static class Gap implements Comparable<Gap> {
		int distance;
		SkeletonMinutia end1;
		SkeletonMinutia end2;
		@Override public int compareTo(Gap other) {
			return Integer.compare(distance, other.distance);
		}
	}
	private void removeGaps() {
		PriorityQueue<Gap> queue = new PriorityQueue<>();
		// 모든 미누셔에 대해
		for (SkeletonMinutia end1 : minutiae)
			// 리지 개수가 1이고 리지의 길이가 일정수준 이상일 때
			if (end1.ridges.size() == 1 && end1.ridges.get(0).points.size() >= Parameters.shortestJoinedEnding)
				// 모든 미누셔에 대해
				for (SkeletonMinutia end2 : minutiae)
					// 다른 미누셔이고 길이가 일정수준 이상이고
					// 갭리미트 보다 작으면 갭으로 추가한다.
					if (end2 != end1 && end2.ridges.size() == 1 && end1.ridges.get(0).end() != end2
						&& end2.ridges.get(0).points.size() >= Parameters.shortestJoinedEnding && isWithinGapLimits(end1, end2)) {
						Gap gap = new Gap();
						gap.distance = end1.position.minus(end2.position).lengthSq();
						gap.end1 = end1;
						gap.end2 = end2;
						queue.add(gap);
					}
		BooleanMap shadow = shadow();
		// 갭이 존재하면
		while (!queue.isEmpty()) {
			// 갭하나를 선택하여
			Gap gap = queue.remove();
			// 갭에 포함된 두 리지가 모두 엔딩이면
			if (gap.end1.ridges.size() == 1 && gap.end2.ridges.size() == 1) {
				Cell[] line = gap.end1.position.lineTo(gap.end2.position);
				if (!isRidgeOverlapping(line, shadow))
					// 두개의 갭을 연결한다.
					addGapRidge(shadow, gap, line);
			}
		}
		removeKnots();
		logger.logRemovedGaps(this);
	}
	private boolean isWithinGapLimits(SkeletonMinutia end1, SkeletonMinutia end2) {
		int distanceSq = end1.position.minus(end2.position).lengthSq();
		// 두개의 미누셔가 가까우면
		if (distanceSq <= Integers.sq(Parameters.maxRuptureSize))
			return true;
		// 충분히 멀면 
		if (distanceSq > Integers.sq(Parameters.maxGapSize))
			return false;
		// 갭간의 각도를 계산하고
		double gapDirection = Angle.atan(end1.position, end2.position);
		// 리지를 따라가 끝점이나 22번째 점까지의 각도를 계산
		double direction1 = Angle.atan(end1.position, angleSampleForGapRemoval(end1));
		// 45도보다 크면 리젝트
		if (Angle.distance(direction1, Angle.opposite(gapDirection)) > Parameters.maxGapAngle)
			return false;
		double direction2 = Angle.atan(end2.position, angleSampleForGapRemoval(end2));
		// 45도보다 크면 리젝트
		if (Angle.distance(direction2, gapDirection) > Parameters.maxGapAngle)
			return false;
		return true;
	}
	private Cell angleSampleForGapRemoval(SkeletonMinutia minutia) {
		SkeletonRidge ridge = minutia.ridges.get(0);
		// 22번째 리지 점 또는 마지막 점 중 가까운 것을 선택
		if (Parameters.gapAngleOffset < ridge.points.size())
			return ridge.points.get(Parameters.gapAngleOffset);
		else
			return ridge.end().position;
	}
	private boolean isRidgeOverlapping(Cell[] line, BooleanMap shadow) {
		for (int i = Parameters.toleratedGapOverlap; i < line.length - Parameters.toleratedGapOverlap; ++i)
			if (shadow.get(line[i]))
				return true;
		return false;
	}
	private static void addGapRidge(BooleanMap shadow, Gap gap, Cell[] line) {
		SkeletonRidge ridge = new SkeletonRidge();
		for (Cell point : line)
			ridge.points.add(point);
		ridge.start(gap.end1);
		ridge.end(gap.end2);
		for (Cell point : line)
			shadow.set(point, true);
	}
	private void removeTails() {
		for (SkeletonMinutia minutia : minutiae) {
			if (minutia.ridges.size() == 1 && minutia.ridges.get(0).end().ridges.size() >= 3)
				if (minutia.ridges.get(0).points.size() < Parameters.minTailLength)
					minutia.ridges.get(0).detach();
		}
		removeDots();
		removeKnots();
		logger.logRemovedTails(this);
	}
	private void removeFragments() {
		for (SkeletonMinutia minutia : minutiae)
			if (minutia.ridges.size() == 1) {
				SkeletonRidge ridge = minutia.ridges.get(0);
				if (ridge.end().ridges.size() == 1 && ridge.points.size() < Parameters.minFragmentLength)
					ridge.detach();
			}
		removeDots();
		logger.logRemovedFragments(this);
	}
	private void removeKnots() {
		for (SkeletonMinutia minutia : minutiae) {
			if (minutia.ridges.size() == 2 && minutia.ridges.get(0).reversed != minutia.ridges.get(1)) {
				SkeletonRidge extended = minutia.ridges.get(0).reversed;
				SkeletonRidge removed = minutia.ridges.get(1);
				if (extended.points.size() < removed.points.size()) {
					SkeletonRidge tmp = extended;
					extended = removed;
					removed = tmp;
					extended = extended.reversed;
					removed = removed.reversed;
				}
				extended.points.remove(extended.points.size() - 1);
				for (Cell point : removed.points)
					extended.points.add(point);
				extended.end(removed.end());
				removed.detach();
			}
		}
		removeDots();
	}
	private void addMinutia(SkeletonMinutia minutia) {
		minutiae.add(minutia);
	}
	private void removeMinutia(SkeletonMinutia minutia) {
		minutiae.remove(minutia);
	}
	private BooleanMap shadow() {
		BooleanMap shadow = new BooleanMap(size);
		for (SkeletonMinutia minutia : minutiae) {
			shadow.set(minutia.position, true);
			for (SkeletonRidge ridge : minutia.ridges)
				if (ridge.start().position.y <= ridge.end().position.y)
					for (Cell point : ridge.points)
						shadow.set(point, true);
		}
		return shadow;
	}
	ByteBuffer serialize() {
		ByteBuffer buffer = ByteBuffer.allocate(minutiae.stream().mapToInt(m -> m.serializedSize()).sum());
		for (SkeletonMinutia minutia : minutiae)
			minutia.write(buffer);
		buffer.flip();
		return buffer;
	}
}
