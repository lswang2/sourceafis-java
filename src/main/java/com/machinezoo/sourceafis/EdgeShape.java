// Part of SourceAFIS: https://sourceafis.machinezoo.com
package com.machinezoo.sourceafis;

class EdgeShape {
	private static final int polarCacheBits = 8;
	private static final int polarCacheRadius = 1 << polarCacheBits;
	private static final int[] polarDistance = new int[Integers.sq(polarCacheRadius)];
	private static final double[] polarAngle = new double[Integers.sq(polarCacheRadius)];
	final int length;
	final double referenceAngle;
	final double neighborAngle;
	static {
		// 256x256테이블로 거리를 미리 계산해 놓는다
		// 256x256테이블로 벡터의 방향을 미리 계산해 놓는다
		// 벡터의 방향은 1사분면만 계산하고 나머지는 회전하여 사용한다.
		for (int y = 0; y < polarCacheRadius; ++y)
			for (int x = 0; x < polarCacheRadius; ++x) {
				polarDistance[y * polarCacheRadius + x] = (int)Math.round(Math.sqrt(Integers.sq(x) + Integers.sq(y)));
				if (y > 0 || x > 0)
					polarAngle[y * polarCacheRadius + x] = Angle.atan(new Point(x, y));
				else
					polarAngle[y * polarCacheRadius + x] = 0;
			}
	}
	EdgeShape(int length, double referenceAngle, double neighborAngle) {
		this.length = length;
		this.referenceAngle = referenceAngle;
		this.neighborAngle = neighborAngle;
	}
	EdgeShape(Minutia reference, Minutia neighbor) {
		// 두개의 미누셔간의 거리를 계산
		// 두개의 미누셔의 방향을 에지벡터로 보정함??
		// 이렇게 각도를 보정하면 지문의 회전속성을 완전히 제거 가능해보임

		// neighbor -> reference로의 벡터
		Cell vector = neighbor.position.minus(reference.position);
		// vector는 무조건 1사분면에 위치시키며
		// quadrant를 통해 실제 어느 사분면에 있어야 되는지 표시한다.
		double quadrant = 0;
		int x = vector.x;
		int y = vector.y;
		if (y < 0) {
			// reference가 아래에 있으면
			// 벡터의 방향을 바꿈
			x = -x;
			y = -y;
			quadrant = Math.PI;
		}
		if (x < 0) {
			// reference가 오른쪽에 있으면
			// 벡터 방향을 -90도 돌림
			int tmp = -x;
			x = y;
			y = tmp;
			quadrant += Angle.halfPI;
		}
		// 벡터의 숫자크기를 알아내서 8(polarCacheBits)을 빼둔다
		int shift = 32 - Integer.numberOfLeadingZeros((x | y) >>> polarCacheBits);
		// x,y를 shift만큼하면 데이터는 무조건 1~255사이의 값이 된다.
		// offset은 16bit값
		int offset = (y >> shift) * polarCacheRadius + (x >> shift);
		// 미리 계산된 테이블에서 계산값을 긁어와서 원래 스케일로 복원한다.
		// 약간의 오차가 있지만 특별히 문제가 될 정도는 아님
		length = polarDistance[offset] << shift;
		// 각도는 테이블에서 읽어서 각 사분면으로 회전하여 사용
		double angle = polarAngle[offset] + quadrant;

		referenceAngle = Angle.difference(reference.direction, angle);
		neighborAngle = Angle.difference(neighbor.direction, Angle.opposite(angle));
	}
}
