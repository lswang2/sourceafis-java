// Part of SourceAFIS: https://sourceafis.machinezoo.com
package com.machinezoo.sourceafis;

class Score {
	int matchedMinutiae;
	double matchedMinutiaeScore;
	double matchedFractionOfProbeMinutiae;
	double matchedFractionOfCandidateMinutiae;
	double matchedFractionOfAllMinutiaeScore;
	int matchedEdges;
	double matchedEdgesScore;
	int minutiaeWithSeveralEdges;
	double minutiaeWithSeveralEdgesScore;
	int correctMinutiaTypeCount;
	double correctMinutiaTypeScore;
	double accurateEdgeLengthScore;
	double accurateMinutiaAngleScore;
	double totalScore;
	double shapedScore;
	void compute(MatchBuffer match) {
		// 현재 root로부터 뻗어나가는 모든 minutia의 개수
		matchedMinutiae = match.count;
		// 매치된 minutia개수를 스코어에 반영 
		matchedMinutiaeScore = Parameters.pairCountScore * matchedMinutiae;

		// probe전체 minutia중 매치된 minutia의 비율
		matchedFractionOfProbeMinutiae = match.count / (double)match.probe.minutiae.length;
		// candidate전체 minutia중 매치된 것의 비율
		matchedFractionOfCandidateMinutiae = match.count / (double)match.candidate.minutiae.length;
		// 비율을 이용한 스코어
		matchedFractionOfAllMinutiaeScore = Parameters.pairFractionScore * (matchedFractionOfProbeMinutiae + matchedFractionOfCandidateMinutiae) / 2;

		// supporingEdge를 반영하여 총 매치된 edge를 카운트
		matchedEdges = match.count;
		minutiaeWithSeveralEdges = 0;
		correctMinutiaTypeCount = 0;
		for (int i = 0; i < match.count; ++i) {
			MinutiaPair pair = match.tree[i];
			matchedEdges += pair.supportingEdges;
			// 중복Edge를 추가로 고려
			if (pair.supportingEdges >= Parameters.minSupportingEdges)
				++minutiaeWithSeveralEdges;
			//매치된 minutia의 타잎이 같으면
			if (match.probe.minutiae[pair.probe].type == match.candidate.minutiae[pair.candidate].type)
				++correctMinutiaTypeCount;
		}
		// 총 edge수를 스코어에 반영 
		matchedEdgesScore = Parameters.edgeCountScore * matchedEdges;
		// 중복edge수를 스코어에 반영
		minutiaeWithSeveralEdgesScore = Parameters.supportedCountScore * minutiaeWithSeveralEdges;
		// type일치 개수를 스코어에 반영
		correctMinutiaTypeScore = Parameters.correctTypeScore * correctMinutiaTypeCount;

		//
		int innerDistanceRadius = (int)Math.round(Parameters.distanceErrorFlatness * Parameters.maxDistanceError);
		int innerAngleRadius = (int)Math.round(Parameters.angleErrorFlatness * Parameters.maxAngleError);
		int distanceErrorSum = 0;
		int angleErrorSum = 0;
		for (int i = 1; i < match.count; ++i) {
			MinutiaPair pair = match.tree[i];
			// root방향의 edge에대하여
			EdgeShape probeEdge = new EdgeShape(match.probe.minutiae[pair.probeRef], match.probe.minutiae[pair.probe]);
			EdgeShape candidateEdge = new EdgeShape(match.candidate.minutiae[pair.candidateRef], match.candidate.minutiae[pair.candidate]);
			//오차가 innerDistanceRadius보다 작으면 innerDistanceRadius를 추가
			distanceErrorSum += Math.max(innerDistanceRadius, Math.abs(probeEdge.length - candidateEdge.length));
			// angle오차도 각각 추가
			angleErrorSum += Math.max(innerAngleRadius, Angle.distance(probeEdge.referenceAngle, candidateEdge.referenceAngle));
			angleErrorSum += Math.max(innerAngleRadius, Angle.distance(probeEdge.neighborAngle, candidateEdge.neighborAngle));
		}
		accurateEdgeLengthScore = 0;
		accurateMinutiaAngleScore = 0;
		if (match.count >= 2) {
			// 엣지들의 거리 오차를 총가능오차에서의 비율을 계산
			double pairedDistanceError = Parameters.maxDistanceError * (match.count - 1);
			accurateEdgeLengthScore = Parameters.distanceAccuracyScore * (pairedDistanceError - distanceErrorSum) / pairedDistanceError;
			// 엣지들의 각도 오차를 총가능오차에서의 비율을 계산
			double pairedAngleError = Parameters.maxAngleError * (match.count - 1) * 2;
			accurateMinutiaAngleScore = Parameters.angleAccuracyScore * (pairedAngleError - angleErrorSum) / pairedAngleError;
		}
		totalScore = matchedMinutiaeScore
			+ matchedFractionOfAllMinutiaeScore
			+ minutiaeWithSeveralEdgesScore
			+ matchedEdgesScore
			+ correctMinutiaTypeScore
			+ accurateEdgeLengthScore
			+ accurateMinutiaAngleScore;
		// 스코어를 쉐이핑함
		shapedScore = shape(totalScore);
	}
	private static double shape(double raw) {
		if (raw < Parameters.thresholdMaxFMR)
			return 0;
		if (raw < Parameters.thresholdFMR2)
			return interpolate(raw, Parameters.thresholdMaxFMR, Parameters.thresholdFMR2, 0, 3);
		if (raw < Parameters.thresholdFMR10)
			return interpolate(raw, Parameters.thresholdFMR2, Parameters.thresholdFMR10, 3, 7);
		if (raw < Parameters.thresholdFMR100)
			return interpolate(raw, Parameters.thresholdFMR10, Parameters.thresholdFMR100, 10, 10);
		if (raw < Parameters.thresholdFMR1000)
			return interpolate(raw, Parameters.thresholdFMR100, Parameters.thresholdFMR1000, 20, 10);
		if (raw < Parameters.thresholdFMR10_000)
			return interpolate(raw, Parameters.thresholdFMR1000, Parameters.thresholdFMR10_000, 30, 10);
		if (raw < Parameters.thresholdFMR100_000)
			return interpolate(raw, Parameters.thresholdFMR10_000, Parameters.thresholdFMR100_000, 40, 10);
		return (raw - Parameters.thresholdFMR100_000) / (Parameters.thresholdFMR100_000 - Parameters.thresholdFMR100) * 30 + 50;
	}
	private static double interpolate(double raw, double min, double max, double start, double length) {
		return (raw - min) / (max - min) * length + start;
	}
}
