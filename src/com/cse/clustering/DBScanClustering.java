package com.cse.clustering;

import static com.cse.clustering.ClusteringAlgorithms.*;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.TreeMap;

public class DBScanClustering {
	public static void implementDBScan(double eps, int minPts) {

		for (Integer i : genes.keySet()) {
			if (visitedMap.containsKey(i)) {
				if (!visitedMap.get(i)) {
					ArrayList<Integer> cluster = new ArrayList<>();
					visitedMap.put(i, true);
					ArrayList<Integer> neighbors = regionQuery(i, eps);
					if (neighbors.size() >= minPts) {
						expandCluster(i, neighbors, cluster, eps, minPts);
						DBScanClusters.add(cluster);
					}
				}
			}
		}
	}

	public static ArrayList<Integer> regionQuery(int geneID, double eps) {
		ArrayList<Integer> list = new ArrayList<>();
		ArrayList<Double> expOneValues = genes.get(geneID);
		for (Integer i : genes.keySet()) {
			ArrayList<Double> expTwoValues = genes.get(i);
			double distance = calculateEuclideanDistance(expOneValues, expTwoValues);
			if (distance <= eps) {
				list.add(i);
			}
		}
		list.add(geneID);
		return list;
	}

	public static void expandCluster(int geneId, ArrayList<Integer> neighbors, ArrayList<Integer> cluster, double eps,
			int minPts) {
		cluster.add(geneId);
		genesInsertedIntocluster.add(geneId);
		for (int i = 0; i < neighbors.size(); i++) {
			if (!visitedMap.get(neighbors.get(i))) {
				visitedMap.put(neighbors.get(i), true);
				ArrayList<Integer> newNeigh = regionQuery(neighbors.get(i), eps);
				if (newNeigh.size() >= minPts) {
					neighbors.addAll(newNeigh);
				}
			}

			if (!genesInsertedIntocluster.contains(neighbors.get(i))) {
				cluster.add(neighbors.get(i));
				genesInsertedIntocluster.add(neighbors.get(i));
			}
		}
	}

	public static void fillVisitedMap() {
		for (Integer i : genes.keySet()) {
			visitedMap.put(i, false);
		}
	}

	public static void generateDBScanMatrix() {
		dBScanMatrix = new int[groundTruthMatrix.length][groundTruthMatrix.length];
		for (ArrayList<Integer> list : DBScanClusters) {
			for (int m = 0; m < list.size() - 1; m++) {
				for (int n = m + 1; n < list.size(); n++) {
					dBScanMatrix[list.get(m)][list.get(n)] = 1;
					dBScanMatrix[list.get(n)][list.get(m)] = 1;
				}
			}
		}
	}

	public static double calculateIndicesDBScan() {
		generateDBScanMatrix();
		int m11 = 0;
		int m10 = 0;
		int m01 = 0;
		int m00 = 0;

		for (int i = 1; i < groundTruthMatrix.length; i++) {
			for (int j = 1; j < groundTruthMatrix.length; j++) {

				if (groundTruthMatrix[i][j] == 1 && dBScanMatrix[i][j] == 1) {
					m11++;
				} else if (groundTruthMatrix[i][j] == 0 && dBScanMatrix[i][j] == 0) {
					m00++;
				} else if (groundTruthMatrix[i][j] == 0 && dBScanMatrix[i][j] == 1) {
					m01++;
				} else if (groundTruthMatrix[i][j] == 1 && dBScanMatrix[i][j] == 0) {
					m10++;
				}
			}
		}

		System.out.println("Jaccard Index: " + m11 / ((m11 + m10 + m01) * 1.0));
		System.out.println("Rand Index: " + (m11 + m00) / ((m00 + m11 + m10 + m01) * 1.0));
		return 0;
	}

	public static void writeToFileAssignedList(ArrayList<ArrayList<Integer>> DBScanClusters, String fileName) {
		try {
			TreeMap<Integer, Integer> map = new TreeMap<>();
			for (int i = 1; i <= genes.size(); i++) {
				map.put(i, 0);
			}

			PrintWriter writer = new PrintWriter("Data/" + fileName + ".txt", "UTF-8");
			for (int cluster = 0; cluster < DBScanClusters.size(); cluster++) {
				for (int i : DBScanClusters.get(cluster)) {
					map.put(i, cluster + 1);
				}
			}

			for (int i : map.values()) {
				writer.println(i);
			}

			writer.close();
		} catch (Exception e) {
			// do something
		}
	}

	public static double calculateSilhouetteCoeff() {
		double silhouette = 0;
		for (Integer geneId : genes.keySet()) {
			int clusterID = getClusterID(geneId);
			if (clusterID != -1) {
				double a = calculateAvgInternalDistanceForSilhouette(geneId, clusterID);
				double b = calculateMinExternalDistanceForSilhouette(geneId, clusterID);
				if (a >= b) {
					silhouette += (b / a) - 1;
				} else {
					silhouette += 1 - (a / b);
				}
			}
		}
		silhouette = silhouette / genes.size();
		System.out.println("Correlation Index (Internal): " + silhouette);
		return silhouette;
	}

	private static int getClusterID(int geneID) {
		for (int i = 0; i < DBScanClusters.size(); i++) {
			if (DBScanClusters.get(i).contains(geneID)) {
				return i;
			}
		}
		return -1;
	}

	private static double calculateAvgInternalDistanceForSilhouette(int geneID, int clusterId) {
		double sum = 0.0;
		for (int g : DBScanClusters.get(clusterId)) {
			sum += calculateEuclideanDistance(genes.get(geneID), genes.get(g));
		}
		return sum / DBScanClusters.get(clusterId).size();
	}

	private static double calculateMinExternalDistanceForSilhouette(int geneID, int clusterId) {
		double min = Double.MAX_VALUE;
		for (int i = 0; i < DBScanClusters.size(); i++) {
			if (i != clusterId) {
				double sum = 0.0;
				for (int g : DBScanClusters.get(i)) {
					sum += calculateEuclideanDistance(genes.get(geneID), genes.get(g));

				}
				sum = sum / DBScanClusters.get(i).size();
				min = Math.min(sum, min);
			}

		}
		return min;
	}

}
