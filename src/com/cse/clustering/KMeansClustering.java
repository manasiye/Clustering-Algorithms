package com.cse.clustering;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.TreeMap;

import static com.cse.clustering.ClusteringAlgorithms.*;

public class KMeansClustering {

	public static void implementKMeans(int k, int numOfIterations, int[] array) {
		if (array.length == 0) {
			initializeCentroid(k);
		} else {
			// System.out.println("Nothing Random");
			initializeWithCentroids(array, k);
		}
		while (numOfIterations-- > 0) {
			prevCluster = new ArrayList<>(clusters);
			clusters.clear();
			for (int i = 0; i < k; i++)
				clusters.add(new ArrayList<Integer>());
			for (Integer g : genes.keySet()) {
				ArrayList<Double> list = new ArrayList<>(genes.get(g));
				int clusterVal = -1;
				Double result = Double.MAX_VALUE;
				for (Integer cent : centroid.keySet()) {
					Double r = calculateEuclideanDistance(centroid.get(cent), list);
					if (r < result) {
						result = r;
						clusterVal = cent;
					}
				}
				clusters.get(clusterVal).add(g);
			}
			// System.out.println(numOfIterations);
			if (AreTwoClustersEqual()) {
				return;
			}
			refillCentroidsWithNewMeans();
		}
	}

	public static void initializeWithCentroids(int[] array, int k) {
		for (int i = 0; i < k; i++) {
			centroid.put(i, genes.get(array[i]));
		}
	}

	public static void initializeCentroid(int k) {
		int size = genes.size();
		Random randomNumberGenerator = new Random();
		for (int i = 0; i < k; i++) {
			int randomNumber = randomNumberGenerator.nextInt(size);
			if (randomNumber == 0) {
				randomNumber = randomNumberGenerator.nextInt(size);
			}
			centroid.put(i, genes.get(randomNumber));
		}
	}

	private static void refillCentroidsWithNewMeans() {
		centroid.clear();
		int numberOfCentroid = clusters.size();
		for (int i = 0; i < numberOfCentroid; i++) {
			ArrayList<Integer> listOfGeneIds = clusters.get(i);
			int expLength = genes.get(1/* listOfGeneIds.get(0) */).size();
			int numberOfGenesInCluster = listOfGeneIds.size();
			ArrayList<Double> listOfExpValues = new ArrayList<>();
			for (int j = 0; j < expLength; j++) {
				double sum = 0;
				for (int k = 0; k < numberOfGenesInCluster; k++) {
					sum = sum + genes.get(listOfGeneIds.get(k)).get(j);
				}
				sum = sum / numberOfGenesInCluster;
				listOfExpValues.add(sum);
			}
			centroid.put(i, listOfExpValues);
		}

	}

	private static boolean AreTwoClustersEqual() {
		if (prevCluster.size() == 0 || clusters.size() == 0)
			return false;
		for (int i = 0; i < clusters.size(); i++) {
			ArrayList<Integer> list1 = clusters.get(i);
			ArrayList<Integer> list2 = prevCluster.get(i);
			if (!equalLists(list1, list2))
				return false;
		}

		return true;
	}

	public static boolean equalLists(ArrayList<Integer> one, ArrayList<Integer> two) {
		if (one == null || two == null) {
			return false;
		}

		if (one.size() != two.size()) {
			return false;
		}

		one = new ArrayList<Integer>(one);
		two = new ArrayList<Integer>(two);

		Collections.sort(one);
		Collections.sort(two);
		return one.equals(two);
	}

	public static void generateKMeansMatrix() {
		kMeansMatrix = new int[groundTruthMatrix.length][groundTruthMatrix.length];
		for (ArrayList<Integer> list : clusters) {
			for (int m = 0; m < list.size() - 1; m++) {
				for (int n = m + 1; n < list.size(); n++) {
					kMeansMatrix[list.get(m)][list.get(n)] = 1;
					kMeansMatrix[list.get(n)][list.get(m)] = 1;
				}
			}
		}
	}

	public static double calculateIndicesKMeans() {
		int m11 = 0;
		int m10 = 0;
		int m01 = 0;
		int m00 = 0;

		for (int i = 1; i < groundTruthMatrix.length; i++) {
			for (int j = 1; j < groundTruthMatrix.length; j++) {

				if (groundTruthMatrix[i][j] == 1 && kMeansMatrix[i][j] == 1) {
					m11++;
				} else if (groundTruthMatrix[i][j] == 0 && kMeansMatrix[i][j] == 0) {
					m00++;
				} else if (groundTruthMatrix[i][j] == 0 && kMeansMatrix[i][j] == 1) {
					m01++;
				} else if (groundTruthMatrix[i][j] == 1 && kMeansMatrix[i][j] == 0) {
					m10++;
				}
			}
		}

		System.out.println("Jaccard Index: " + m11 / ((m11 + m10 + m01) * 1.0));
		System.out.println("Rand Index: " + (m11 + m00) / ((m00 + m11 + m10 + m01) * 1.0));
		return 0;

	}

	public static double calculateSilhouetteCoeff() {
		double silhouette = 0;
		for (Integer geneId : genes.keySet()) {
			int clusterID = getClusterID(geneId);
			double a = calculateAvgInternalDistanceForSilhouette(geneId, clusterID);
			double b = calculateMinExternalDistanceForSilhouette(geneId, clusterID);
			if (a >= b) {
				silhouette += (b / a) - 1;
			} else {
				silhouette += 1 - (a / b);
			}
		}
		silhouette = silhouette / genes.size();
		System.out.println("Correlation Index (Internal): " + silhouette);
		return silhouette;
	}

	private static int getClusterID(int geneID) {
		for (int i = 0; i < clusters.size(); i++) {
			for (int g : clusters.get(i)) {
				if (g == geneID)
					return i;
			}
		}
		return -1;
	}

	private static double calculateAvgInternalDistanceForSilhouette(int geneID, int clusterId) {
		double sum = 0.0;
		for (int g : clusters.get(clusterId)) {
			sum += calculateEuclideanDistance(genes.get(geneID), genes.get(g));
		}
		return sum / clusters.get(clusterId).size();
	}

	private static double calculateMinExternalDistanceForSilhouette(int geneID, int clusterId) {
		double min = Double.MAX_VALUE;
		for (int i = 0; i < clusters.size(); i++) {
			if (i != clusterId) {
				double sum = 0.0;
				for (int g : clusters.get(i)) {
					sum += calculateEuclideanDistance(genes.get(geneID), genes.get(g));

				}
				sum = sum / clusters.get(i).size();
				min = Math.min(sum, min);
			}

		}
		return min;
	}

	public static void writeToFileClusterList() {
		try {
			PrintWriter writer = new PrintWriter("/home/prasad/kMeansResults.txt", "UTF-8");
			for (ArrayList<Integer> list : clusters) {
				for (int i : list) {
					writer.print(i + " ");
				}
				writer.println();
			}

			writer.close();
		} catch (Exception e) {
			// do something
		}
	}

	public static void writeToFileAssignedList(ArrayList<ArrayList<Integer>> clusters, String fileName) {
		try {
			TreeMap<Integer, Integer> map = new TreeMap<>();

			PrintWriter writer = new PrintWriter("Data/" + fileName + ".txt", "UTF-8");
			for (int cluster = 0; cluster < clusters.size(); cluster++) {
				for (int i : clusters.get(cluster)) {
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
}
