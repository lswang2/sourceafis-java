package sourceafis.collections;

import sourceafis.scalars.*;

public class Histogram {
	public final int width;
	public final int height;
	public final int depth;
	public final int[] array;
	public Histogram(int width, int height, int depth) {
		this.width = width;
		this.height = height;
		this.depth = depth;
		array = new int[width * height * depth];
	}
	public int constrain(int z) {
		return Math.max(0, Math.min(depth - 1, z));
	}
	public int get(int x, int y, int z) {
		return array[(y * width + x) * depth + z];
	}
	public int get(Cell at, int z) {
		return get(at.x, at.y, z);
	}
	public void set(int x, int y, int z, int value) {
		array[(y * width + x) * depth + z] = value;
	}
	public void set(Cell at, int z, int value) {
		set(at.x, at.y, z, value);
	}
	public void add(int x, int y, int z, int value) {
		array[(y * width + x) * depth + z] += value;
	}
	public void add(Cell at, int z, int value) {
		add(at.x, at.y, z, value);
	}
	public void increment(int x, int y, int z) {
		add(x, y, z, 1);
	}
	public void increment(Cell at, int z) {
		increment(at.x, at.y, z);
	}
}