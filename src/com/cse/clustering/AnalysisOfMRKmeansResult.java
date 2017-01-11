package com.cse.clustering;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

public class AnalysisOfMRKmeansResult {
	public static HashMap<Integer, ArrayList<Double>> genes = new HashMap<>();
	public static HashMap<Integer, Integer> groundTruth = new HashMap<>();
	public static TreeMap<Integer, Integer> mapGeneIDClusterID = new TreeMap<Integer, Integer>();
	public static int[][] groundTruthMatrix;
	public static int[][] MRKmeansMatrix;

	public AnalysisOfMRKmeansResult(String fileName) {

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
				if ((groundTruth.get(i) == groundTruth.get(j)) && (groundTruth.get(i) != -1)) {
					groundTruthMatrix[i][j] = 1;
				} else {
					groundTruthMatrix[i][j] = 0;
				}
			}
		}
	}

	private static void fillMRKmeansMatrix() {
		MRKmeansMatrix = new int[groundTruthMatrix.length][groundTruthMatrix.length];
		for (int i : mapGeneIDClusterID.keySet()) {
			for (int j : mapGeneIDClusterID.keySet()) {
				if (mapGeneIDClusterID.get(i) == mapGeneIDClusterID.get(j)) {
					MRKmeansMatrix[i][j] = 1;
				}
			}
		}
	}

	private static void writeToFileMROutput(String fileName) {
		String line = null;
		try {
			FileReader fileReader = new FileReader(fileName);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			int clusterID = 0;
			while ((line = bufferedReader.readLine()) != null) {
				String[] lineSplit = line.split("\\[")[1].split(",");
				int len = lineSplit.length;
				for (int i = 0; i < len; i++) {
					if (i == len - 1) {
						mapGeneIDClusterID.put(Integer.parseInt(lineSplit[i].replace("]", "").trim()), clusterID);
					} else {
						mapGeneIDClusterID.put(Integer.parseInt(lineSplit[i].trim()), clusterID);
					}
				}
				clusterID++;
			}
			writeToFileForPCA(mapGeneIDClusterID);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void writeToFileForPCA(TreeMap<Integer, Integer> map) {
		String fileName = "mrAnalysis.txt";
		try {
			PrintWriter writer = new PrintWriter("Data/" + fileName, "UTF-8");

			for (int i : map.values()) {
				writer.println(i);
			}

			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static double calculateJaccardRandMR() {

		int m11 = 0;
		int m10 = 0;
		int m01 = 0;
		int m00 = 0;

		for (int i = 1; i < groundTruthMatrix.length; i++) {
			for (int j = 1; j < groundTruthMatrix.length; j++) {

				if (groundTruthMatrix[i][j] == 1 && MRKmeansMatrix[i][j] == 1) {
					m11++;
				} else if (groundTruthMatrix[i][j] == 0 && MRKmeansMatrix[i][j] == 0) {
					m00++;
				} else if (groundTruthMatrix[i][j] == 0 && MRKmeansMatrix[i][j] == 1) {
					m01++;
				} else if (groundTruthMatrix[i][j] == 1 && MRKmeansMatrix[i][j] == 0) {
					m10++;
				}
			}
		}

		System.out.println("Jaccard Index: " + m11 / ((m11 + m10 + m01) * 1.0));
		System.out.println("Rand Index: " + (m11 + m00) / ((m00 + m11 + m10 + m01) * 1.0));
		return 0;

	}

	public static void main(String[] args) {
		AnalysisOfMRKmeansResult analysisOfMRKmeansResult = new AnalysisOfMRKmeansResult("Data/iyer.txt");
		writeToFileMROutput("Data/fileInput.txt");
		// System.out.println(mapGeneIDClusterID);
		fillMRKmeansMatrix();
		calculateJaccardRandMR();
	}

}
