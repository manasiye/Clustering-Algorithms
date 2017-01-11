package com.cse.clustering;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

public class ClusteringAlgorithms {
	public static HashMap<Integer, ArrayList<Double>> genes = new HashMap<>();
	public static HashMap<Integer, Integer> groundTruth = new HashMap<>();
	public static HashMap<Integer, ArrayList<Double>> centroid = new HashMap<>();
	public static ArrayList<ArrayList<Integer>> clusters = new ArrayList<>();
	public static ArrayList<ArrayList<Integer>> prevCluster = new ArrayList<>();
	public static ArrayList<ArrayList<Integer>> HClusters = new ArrayList<>();
	public static ArrayList<ArrayList<Integer>> DBScanClusters = new ArrayList<>();
	public static HashMap<Integer, Boolean> visitedMap = new HashMap<>();
	public static HashSet<Integer> genesInsertedIntocluster = new HashSet<>();
	public static double[][] matrix;
	public static int[] minDistArray;

	public static ArrayList<HashSet<Integer>> clusterH = new ArrayList<>();
	public static int[][] groundTruthMatrix;
	public static int[][] kMeansMatrix;
	public static int[][] hierarchicalMatrix;
	public static int[][] dBScanMatrix;

	public ClusteringAlgorithms(String fileName) {
		String line = null;
		try {
			FileReader fileReader = new FileReader(fileName);
			@SuppressWarnings("resource")
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			while ((line = bufferedReader.readLine()) != null) {
				String[] lineSplit = line.trim().split("\\t");
				int len = lineSplit.length;
				// System.out.println("0:::" + lineSplit[0]);
				int geneID = Integer.parseInt(lineSplit[0].trim());
				groundTruth.put(geneID, Integer.parseInt(lineSplit[1]));
				ArrayList<Double> list = new ArrayList<>();
				for (int i = 2; i < len; i++) {
					list.add(Double.parseDouble(lineSplit[i]));
				}

				genes.put(geneID, list);
			}
			// System.out.println("S" + groundTruth.size());
			generateIncidenceMatrixGroundTruth();

		} catch (FileNotFoundException e) {
			System.out.println("Unable to open file" + fileName);
		} catch (IOException e) {
			System.out.println("Error reading from file" + fileName);
		}
	}

	public static void generateIncidenceMatrixGroundTruth() {
		groundTruthMatrix = new int[groundTruth.size() + 1][groundTruth.size() + 1];
		for (int i = 1; i < groundTruthMatrix.length; i++) {
			for (int j = 1; j < groundTruthMatrix.length; j++) {
				if (groundTruth.get(i) == groundTruth.get(j) && groundTruth.get(i) != -1) {
					groundTruthMatrix[i][j] = 1;
				} else {
					groundTruthMatrix[i][j] = 0;
				}
			}
		}
	}

	public static double calculateEuclideanDistance(ArrayList<Double> list1, ArrayList<Double> list2) {
		int size = list1.size();
		double euclideanDist = 0;
		for (int i = 0; i < size; i++) {
			double val1 = list1.get(i);
			double val2 = list2.get(i);
			double diff = val1 - val2;
			diff = diff * diff;
			euclideanDist += diff;
		}
		return Math.sqrt(euclideanDist);
	}

	public static void main(String[] args) {

		String s1 = "Data/new_dataset_1.txt";
		String s2 = "Data/new_dataset_2.txt";
		String s3 = "Data/cho.txt";
		String s4 = "Data/iyer.txt";

		// <<<<<<<<<<<<<CONFIG>>>>>>>>>>>>>>>>>>>

		ClusteringAlgorithms kmeans = new ClusteringAlgorithms(s4);

		int kNum = 10;
		int numberOfIterations = Integer.MAX_VALUE;
		int[] beginnerArray = {};

		int numOfClustersForHierarchical = 9;

		double eps = 0.8;
		int minPts = 5;

		// <<<<<<<<<<<<<<<<<<<Kmeans Clustering>>>>>>>>>>>>>>>//

		KMeansClustering.implementKMeans(kNum, numberOfIterations, beginnerArray);
		for (ArrayList<Integer> arr : clusters) {
			System.out.println(arr);
		}
		KMeansClustering.generateKMeansMatrix();
		KMeansClustering.calculateIndicesKMeans();
		KMeansClustering.calculateSilhouetteCoeff();
		KMeansClustering.writeToFileAssignedList(clusters, "kMeans");
		System.out.println(centroid);
		// <<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>//

		// <<<<<<<<<<<<<<<<<<<<<<Hierarchical>>>>>>>>>>>>>>>>>>>>>>//

		HierarchicalClustering.fillHCloud(genes.size());
		HierarchicalClustering.fillHClustersMatrix();

		HierarchicalClustering.implementHierarchical(numOfClustersForHierarchical);
		for (HashSet<Integer> s : clusterH) {
			System.out.println(s);
		}
		HierarchicalClustering.generateHierarchicalMatrix();
		HierarchicalClustering.calculateIndicesHierarchical();
		HierarchicalClustering.calculateSilhouetteCoeff();
		HierarchicalClustering.writeToFileAssignedList(clusterH, "Hierarchical");

		// // <<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>//
		//
		// // <<<<<<<<<<<<<<<<<<<<<<DBScan>>>>>>>>>>>>>>>>>>>>>>//
		//
		DBScanClustering.fillVisitedMap();
		DBScanClustering.implementDBScan(eps, minPts);
		System.out.println("Size: " + DBScanClusters.size());
		for (ArrayList<Integer> list : DBScanClusters) {
			System.out.println(list);
		}
		DBScanClustering.calculateIndicesDBScan();
		DBScanClustering.calculateSilhouetteCoeff();
		DBScanClustering.writeToFileAssignedList(DBScanClusters, "DBScan");

	}
}
