// Part of SourceAFIS: https://sourceafis.machinezoo.com
package com.machinezoo.sourceafis;

import java.util.*;
import gnu.trove.map.hash.*;
import gnu.trove.set.hash.*;

class MatchBuffer {
	private static final ThreadLocal<MatchBuffer> local = ThreadLocal.withInitial(MatchBuffer::new);
	FingerprintTransparency transparency = FingerprintTransparency.none;
	ImmutableTemplate probe;
	// integer를 인덱스로 하는 해시맵
	private TIntObjectHashMap<List<IndexedEdge>> edgeHash;
	ImmutableTemplate candidate;
	private MinutiaPair[] pool = new MinutiaPair[1];
	private int pooled;
	private PriorityQueue<MinutiaPair> queue = new PriorityQueue<>(Comparator.comparing(p -> p.distance));
	int count;
	MinutiaPair[] tree;
	private MinutiaPair[] byProbe;
	private MinutiaPair[] byCandidate;
	private MinutiaPair[] roots;
	private final TIntHashSet duplicates = new TIntHashSet();
	private Score score = new Score();
	static MatchBuffer current() {
		return local.get();
	}
	void selectMatcher(ImmutableMatcher matcher) {
		// probe.minutiae는 원본의 minutia의 array
		// probe.edges는 사용하지 않는가?
		// --> tree구성할 때 사용함
		probe = matcher.template;
		if (tree == null || probe.minutiae.length > tree.length) {
			tree = new MinutiaPair[probe.minutiae.length];
			// minutia개수의 자료??
			byProbe = new MinutiaPair[probe.minutiae.length];
		}
		// edgeHash는 에지의 검색을 빠르게 하기 위한 해시맵
		edgeHash = matcher.edgeHash;
	}
	void selectCandidate(ImmutableTemplate template) {
		// candidate.minutiae는 비교할 새 지문의 minutia의 array
		candidate = template;
		if (byCandidate == null || byCandidate.length < candidate.minutiae.length)
			// minutia개수의 자료??
			byCandidate = new MinutiaPair[candidate.minutiae.length];
	}
	double match() {
		try {
			int totalRoots = enumerateRoots();
			transparency.logRootPairs(totalRoots, roots);
			double high = 0;
			int best = -1;
			for (int i = 0; i < totalRoots; ++i) {
				double partial = tryRoot(roots[i]);
				if (partial > high) {
					high = partial;
					best = i;
				}
				clearPairing();
			}
			transparency.logBestMatch(best);
			// 가장 높은 score반환
			return high;
		} catch (Throwable e) {
			local.remove();
			throw e;
		}
	}
	private int enumerateRoots() {
		if (roots == null || roots.length < Parameters.maxTriedRoots)
			roots = new MinutiaPair[Parameters.maxTriedRoots];
		int totalLookups = 0;
		int totalRoots = 0;
		int triedRoots = 0;
		duplicates.clear();
		for (boolean shortEdges : new boolean[] { false, true }) {
			// period : 1               2                       3
			// phase  : 0       1       0       1       2       0
			// refer  : 0 2 ... 1 3 ... 0 3 ... 1 4 ... 2 5 ... 0 4 ...
			for (int period = 1; period < candidate.minutiae.length; ++period) {
				// candidate의 모든 minutia에 대해
				for (int phase = 0; phase <= period; ++phase) {
					// bubble match??
					for (int candidateReference = phase; candidateReference < candidate.minutiae.length; candidateReference += period + 1) {
						
						// overflow??
						int candidateNeighbor = (candidateReference + period) % candidate.minutiae.length;

						EdgeShape candidateEdge = new EdgeShape(candidate.minutiae[candidateReference], candidate.minutiae[candidateNeighbor]);

						// minRootEdgeLength를 기준으로 짧은 것과 긴 것을 구분하여 처리
						if ((candidateEdge.length >= Parameters.minRootEdgeLength) ^ shortEdges) {

							// probe이미지를 이용하여 만든 edgeHash에서 뭔가를 끄집어냄??
							// 이러면 약간의 오차가 있는 경우 검색이 안되지 않을까?
							// --> 해시맵을 구성할 때 오차까지 고려하고 구성하여 해결함
							List<IndexedEdge> matches = edgeHash.get(hashShape(candidateEdge));
							if (matches != null) {
								//비슷한 edge가 있으면
								for (IndexedEdge match : matches) {
									// 매치된 모든 EdgeShape에 대해서
									// 비교하여 오차이내로 일치하면
									if (matchingShapes(match, candidateEdge)) {
										int duplicateKey = (match.reference << 16) | candidateReference;
										// 기존에 검색한 적이 없으면
										// edge가 아니라 minutia에 대해 duplication을 체크하고
										// probe/candidcate의 minutiaPair를 만듬
										if (!duplicates.contains(duplicateKey)) {
											duplicates.add(duplicateKey);
											MinutiaPair pair = allocate();
											pair.probe = match.reference;
											pair.candidate = candidateReference;
											// roots에 추가
											roots[totalRoots] = pair;
											++totalRoots;
										}
										++triedRoots;
										if (triedRoots >= Parameters.maxTriedRoots)
											return totalRoots;
									}
								}
							}
							++totalLookups;
							if (totalLookups >= Parameters.maxRootEdgeLookups)
								return totalRoots;
						}
					}
				}
			}
		}
		return totalRoots;
	}
	private int hashShape(EdgeShape edge) {
		// EdgeShape를 오차를 배제한 값으로 단일 정수로 생성함
		// 양자화오차의 가능성이 남음
		int lengthBin = edge.length / Parameters.maxDistanceError;
		int referenceAngleBin = (int)(edge.referenceAngle / Parameters.maxAngleError);
		int neighborAngleBin = (int)(edge.neighborAngle / Parameters.maxAngleError);
		return (referenceAngleBin << 24) + (neighborAngleBin << 16) + lengthBin;
	}
	private boolean matchingShapes(EdgeShape probe, EdgeShape candidate) {
		// 두 에지가 비슷하면 true
		int lengthDelta = probe.length - candidate.length;
		if (lengthDelta >= -Parameters.maxDistanceError && lengthDelta <= Parameters.maxDistanceError) {
			double complementaryAngleError = Angle.complementary(Parameters.maxAngleError);
			double referenceDelta = Angle.difference(probe.referenceAngle, candidate.referenceAngle);
			if (referenceDelta <= Parameters.maxAngleError || referenceDelta >= complementaryAngleError) {
				double neighborDelta = Angle.difference(probe.neighborAngle, candidate.neighborAngle);
				if (neighborDelta <= Parameters.maxAngleError || neighborDelta >= complementaryAngleError)
					return true;
			}
		}
		return false;
	}
	private double tryRoot(MinutiaPair root) {
		// queue 하나 루트 minutiaPair삽입
		queue.add(root);
		do {
			// tree에 queue의 맨앞 멤버를 추가
			addPair(queue.remove());
			collectEdges();
			skipPaired();
		} while (!queue.isEmpty());
		transparency.logPairing(count, tree);
		// 현재 루트로 검색된 정보를 가지고 score계산
		score.compute(this);
		transparency.logScore(score);
		return score.shapedScore;
	}
	private void clearPairing() {
		for (int i = 0; i < count; ++i) {
			byProbe[tree[i].probe] = null;
			byCandidate[tree[i].candidate] = null;
			release(tree[i]);
			tree[i] = null;
		}
		count = 0;
	}
	private void collectEdges() {
		// 현재 노드의 모든 에지에 대해 검색해서 일치하는 것은
		// tree child로 추가함
		// 이미 있는 것은 supportingEdge++함

		// tree의 바로 직전 멤버를 읽어냄
		MinutiaPair reference = tree[count - 1];
		// probe에서 이웃 edge를 얻어 냄
		NeighborEdge[] probeNeighbors = probe.edges[reference.probe];
		// candicate에서 이웃 edge를 얻어 냄
		NeighborEdge[] candidateNeigbors = candidate.edges[reference.candidate];
		for (MinutiaPair pair : matchPairs(probeNeighbors, candidateNeigbors)) {
			// probe/candidate에서 일치하는 edge들만의 neighbor들을 minutiaPair의 list로 구성하여 처리함 
			// probeRef는 probe에서 해당 기준minutia의 인덱스
			pair.probeRef = reference.probe;
			// candidateRef는 candidatee에서 해당 기준minutia의 인덱스
			pair.candidateRef = reference.candidate;
			// 이전에 addPair 되지 않았던 멤버들이면 queue에 넣고
			if (byCandidate[pair.candidate] == null && byProbe[pair.probe] == null)
				queue.add(pair);
			else {
				// probe의 트리에는 들어 있고
				// 그 상대 candidate가 서로 일치하면
				// 서포팅 엣지를 카운트업 한다.
				if (byProbe[pair.probe] != null && byProbe[pair.probe].candidate == pair.candidate)
					addSupportingEdge(pair);
				// 해당 pair를 해제한다.
				release(pair);
			}
		}
	}
	private List<MinutiaPair> matchPairs(NeighborEdge[] probeStar, NeighborEdge[] candidateStar) {
		double complementaryAngleError = Angle.complementary(Parameters.maxAngleError);
		List<MinutiaPair> results = new ArrayList<>();
		int start = 0;
		int end = 0;
		// 이웃 엣지들이 길이로 정렬되어 있는 것으로 추축됨
		for (int candidateIndex = 0; candidateIndex < candidateStar.length; ++candidateIndex) {
			NeighborEdge candidateEdge = candidateStar[candidateIndex];
			// start,end를 candidate와 거리오차 이내의 것만을 대상으로 한다.
			while (start < probeStar.length && probeStar[start].length < candidateEdge.length - Parameters.maxDistanceError)
				++start;
			if (end < start)
				end = start;
			while (end < probeStar.length && probeStar[end].length <= candidateEdge.length + Parameters.maxDistanceError)
				++end;
			for (int probeIndex = start; probeIndex < end; ++probeIndex) {
				NeighborEdge probeEdge = probeStar[probeIndex];
				double referenceDiff = Angle.difference(probeEdge.referenceAngle, candidateEdge.referenceAngle);
				if (referenceDiff <= Parameters.maxAngleError || referenceDiff >= complementaryAngleError) {
					double neighborDiff = Angle.difference(probeEdge.neighborAngle, candidateEdge.neighborAngle);
					if (neighborDiff <= Parameters.maxAngleError || neighborDiff >= complementaryAngleError) {
						// edge가 일치하는 것들로만 리스트를 구성한다.
						MinutiaPair pair = allocate();
						pair.probe = probeEdge.neighbor;
						pair.candidate = candidateEdge.neighbor;
						// distance는 나중에 queue에서 뽑혀나오는 순서를 결정한다.
						pair.distance = candidateEdge.length;
						results.add(pair);
					}
				}
			}
		}
		return results;
	}
	private void skipPaired() {
		// 현재 큐에서 짧은 에지들은 이미 커버되었을 가능성이 높으므로
		// 미리 검사하여 supporingEdge++로 처리한다.

		// queue가 비어있지 않고
		// queue에 있는 가장 짧은 edge의 네이버가 이미 커버되었던 minutia이면
		while (!queue.isEmpty() && (byProbe[queue.peek().probe] != null || byCandidate[queue.peek().candidate] != null)) {
			// 해당 minutiaPair를 queue에서 제거한 후
			MinutiaPair pair = queue.remove();
			// probe의 트리에는 들어 있고
			// 그 상대 candidate가 서로 일치하면
			// 서포팅 엣지를 카운트업 한다.
			if (byProbe[pair.probe] != null && byProbe[pair.probe].candidate == pair.candidate)
				addSupportingEdge(pair);
			release(pair);
		}
	}
	private void addPair(MinutiaPair pair) {
		// tree 구성
		tree[count] = pair;
		// byProbe는 초기에는 모두 null이었다가 차츰 추가됨
		byProbe[pair.probe] = pair;
		byCandidate[pair.candidate] = pair;
		++count;
	}
	private void addSupportingEdge(MinutiaPair pair) {
		++byProbe[pair.probe].supportingEdges;
		++byProbe[pair.probeRef].supportingEdges;
		transparency.logSupportingEdge(pair);
	}
	private MinutiaPair allocate() {
		if (pooled > 0) {
			--pooled;
			MinutiaPair pair = pool[pooled];
			pool[pooled] = null;
			return pair;
		} else
			return new MinutiaPair();
	}
	private void release(MinutiaPair pair) {
		if (pooled >= pool.length)
			pool = Arrays.copyOf(pool, 2 * pool.length);
		pair.probe = 0;
		pair.candidate = 0;
		pair.probeRef = 0;
		pair.candidateRef = 0;
		pair.distance = 0;
		pair.supportingEdges = 0;
		pool[pooled] = pair;
	}
}
