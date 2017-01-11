package com.cse.clustering;

import static com.cse.clustering.ClusteringAlgorithms.*;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.TreeMap;

public class HierarchicalClustering {

	public static void implementHierarchical(int k) {

		for (int s = 0; s < genes.size(); s++) {
			int num1 = findSmallestInMinArray(1);

			int num2 = minDistArray[num1];
			insertIntoHCloud(num1, num2);
			overwriteRowWithMinimum(num1, num2);

			for (int i = 1; i <= genes.size(); i++) {
				matrix[num2][i] = matrix[i][num2] = Double.POSITIVE_INFINITY;
			}

			updateMinArray(num1, num2);
			if (clusterH.size() <= k) {
				return;
			}
			// System.out.println(i1 + " " + i2);
		}
	}

	private static void overwriteRowWithMinimum(int num1, int num2) {
		for (int j = 1; j <= genes.size(); j++) {
			if (matrix[num2][j] < matrix[num1][j])
				matrix[num1][j] = matrix[j][num1] = matrix[num2][j];
		}
		matrix[num1][num1] = Double.POSITIVE_INFINITY;
	}

	private static void insertIntoHCloud(int val1, int val2) {
		int v1 = -1;
		int v2 = -1;
		for (int i = 0; i < clusterH.size(); i++) {
			if (clusterH.get(i).contains(val1)) {
				v1 = i;
			}
			if (clusterH.get(i).contains(val2)) {
				v2 = i;
			}
		}
		if (v1 != -1 && v2 != -1) {
			// merge to sets
			HashSet<Integer> set = clusterH.get(v2);
			clusterH.get(v1).addAll(set);
			clusterH.remove(v2);

		} else if (v1 == -1 && v2 != -1) {
			clusterH.get(v2).add(val1);
		} else if (v2 == -1 && v1 != -1) {
			clusterH.get(v1).add(val2);
		} else {
			clusterH.add(new HashSet<Integer>());
			clusterH.get(clusterH.size() - 1).add(val1);
			clusterH.get(clusterH.size() - 1).add(val2);
		}
	}

	public static void fillHCloud(int size) {
		for (int i = 1; i <= size; i++) {
			HashSet<Integer> set = new HashSet<Integer>();
			set.add(i);
			clusterH.add(set);
		}
	}

	public static void fillHClustersMatrix() {
		int genesCount = genes.size();
		matrix = new double[genesCount + 1][genesCount + 1];
		minDistArray = new int[genesCount + 1];
		Arrays.fill(minDistArray, 1);

		for (int i = 1; i <= genesCount; i++) {
			for (int j = 1; j <= genesCount; j++) {
				if (i == j) {
					matrix[i][j] = Double.POSITIVE_INFINITY;
				} else {
					ArrayList<Double> list1 = genes.get(i);
					ArrayList<Double> list2 = genes.get(j);
					double dist = calculateEuclideanDistance(list1, list2);
					matrix[i][j] = dist;
				}
				if (matrix[i][j] < matrix[i][minDistArray[i]]) {
					minDistArray[i] = j;
				}
			}
		}
	}

	public static void generateHierarchicalMatrix() {
		hierarchicalMatrix = new int[groundTruthMatrix.length][groundTruthMatrix.length];
		for (HashSet<Integer> list : clusterH) {
			ArrayList<Integer> temp = new ArrayList<>(list);
			for (int m = 0; m < temp.size() - 1; m++) {
				for (int n = m + 1; n < temp.size(); n++) {
					hierarchicalMatrix[temp.get(m)][temp.get(n)] = 1;
					hierarchicalMatrix[temp.get(m)][temp.get(m)] = 1;
				}
			}
		}
	}

	public static double calculateIndicesHierarchical() {
		int m11 = 0;
		int m10 = 0;
		int m01 = 0;
		int m00 = 0;

		for (int i = 1; i < groundTruthMatrix.length; i++) {
			for (int j = 1; j < groundTruthMatrix.length; j++) {

				if (groundTruthMatrix[i][j] == 1 && hierarchicalMatrix[i][j] == 1) {
					m11++;
				} else if (groundTruthMatrix[i][j] == 0 && hierarchicalMatrix[i][j] == 0) {
					m00++;
				} else if (groundTruthMatrix[i][j] == 0 && hierarchicalMatrix[i][j] == 1) {
					m01++;
				} else if (groundTruthMatrix[i][j] == 1 && hierarchicalMatrix[i][j] == 0) {
					m10++;
				}
			}
		}

		System.out.println("Jaccard Index: " + m11 / ((m11 + m10 + m01) * 1.0));
		System.out.println("Rand Index: " + (m11 + m00) / ((m00 + m11 + m10 + m01) * 1.0));
		return 0;

	}

	public static void writeToFileAssignedList(ArrayList<HashSet<Integer>> clusters, String fileName) {
		try {
			TreeMap<Integer, Integer> map = new TreeMap<>();

			PrintWriter writer = new PrintWriter("Data/" + fileName + ".txt", "UTF-8");
			for (int cluster = 0; cluster < clusters.size(); cluster++) {
				for (int i : clusters.get(cluster)) {
					if (i != 0) {
						map.put(i, cluster + 1);
					}
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
			double a = calculateAvgInternalDistanceForSilhouette(geneId, clusterID);
			double b = calculateMinExternalDistanceForSilhouette(geneId, clusterID);
			double tmp = 0.0;
			if (a >= b) {
				tmp += (b / a) - 1;
			} else {
				tmp += 1 - (a / b);
			}

			if (tmp > 0) {
				silhouette += tmp;
			}
		}
		silhouette = silhouette / genes.size();
		System.out.println("Correlation Index (Internal): " + silhouette);
		return silhouette;
	}

	private static int getClusterID(Integer geneId) {
		for (int i = 0; i < clusterH.size(); i++) {
			if (clusterH.get(i).contains(geneId)) {
				return i;
			}
		}
		return -1;
	}

	private static double calculateAvgInternalDistanceForSilhouette(Integer geneId, int clusterID) {
		double sum = 0.0;
		for (Integer i : clusterH.get(clusterID)) {
			sum += calculateEuclideanDistance(genes.get(i), genes.get(geneId));
		}
		return sum / clusterH.get(clusterID).size();
	}

	private static double calculateMinExternalDistanceForSilhouette(Integer geneId, int clusterID) {
		double min = Double.MAX_VALUE;
		for (int i = 0; i < clusterH.size(); i++) {
			if (i != clusterID) {
				double sum = 0.0;
				for (int g : clusterH.get(i)) {
					sum += calculateEuclideanDistance(genes.get(geneId), genes.get(g));

				}
				sum = sum / clusterH.get(i).size();
				min = Math.min(sum, min);
			}
		}
		return min;
	}

	public static int findSmallestInMinArray(int i1) {
		for (int i = 1; i <= genes.size(); i++) {
			if (matrix[i][minDistArray[i]] < matrix[i1][minDistArray[i1]])
				i1 = i;
		}
		return i1;
	}

	public static void updateMinArray(int num1, int num2) {
		for (int j = 1; j <= genes.size(); j++) {
			if (minDistArray[j] == num2)
				minDistArray[j] = num1;
			if (matrix[num1][j] < matrix[num1][minDistArray[num1]])
				minDistArray[num1] = j;
		}
	}

}
