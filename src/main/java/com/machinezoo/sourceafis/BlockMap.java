// Part of SourceAFIS: https://sourceafis.machinezoo.com
package com.machinezoo.sourceafis;

class BlockMap {
	final Cell pixels;
	final BlockGrid primary;
	final BlockGrid secondary;
	BlockMap(int width, int height, int maxBlockSize) {
		pixels = new Cell(width, height);
		// 블록크기를 maxBlockSize로 정하고 
		// 블록 개수를 계산하여 BlockGrid를 생성 
		primary = new BlockGrid(new Cell(
			Integers.roundUpDiv(pixels.x, maxBlockSize),
			Integers.roundUpDiv(pixels.y, maxBlockSize)));
		// BlockGrid의 경계선을 모두 계산한다.
		// 모든 Grid가 동일 크기가 되지는 않는다.
		for (int y = 0; y <= primary.blocks.y; ++y)
			primary.y[y] = y * pixels.y / primary.blocks.y;
		for (int x = 0; x <= primary.blocks.x; ++x)
			primary.x[x] = x * pixels.x / primary.blocks.x;
		// 기존 크기 +1인 그리드
		// primary와 50% overlap되는 grid
		secondary = new BlockGrid(primary.corners);
		secondary.y[0] = 0;
		for (int y = 0; y < primary.blocks.y; ++y)
			secondary.y[y + 1] = primary.block(0, y).center().y;
		secondary.y[secondary.blocks.y] = pixels.y;
		secondary.x[0] = 0;
		for (int x = 0; x < primary.blocks.x; ++x)
			secondary.x[x + 1] = primary.block(x, 0).center().x;
		secondary.x[secondary.blocks.x] = pixels.x;
	}
}
