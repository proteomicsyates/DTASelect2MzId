package edu.scripps.yates.dtaselect2mzid.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

public class MS2Reader {
	private final static Logger log = Logger.getLogger(MS2Reader.class);
	private final File ms2File;
	private Map<String, Integer> spectrumIndexByScanNumberMap;
	private Map<String, Double> rtByScanNumberMap;

	public MS2Reader(File ms2File) {
		this.ms2File = ms2File;
		if (!ms2File.exists()) {
			throw new IllegalArgumentException("File not found at '" + ms2File.getAbsolutePath() + "'");
		}
	}

	public Integer getSpectrumIndexByScan(String scanNumber) {
		process();
		if (spectrumIndexByScanNumberMap.containsKey(scanNumber)) {
			return spectrumIndexByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	private void process() {
		if (spectrumIndexByScanNumberMap == null || rtByScanNumberMap == null) {
			log.info("Reading MS2 file: " + ms2File.getAbsolutePath());
			spectrumIndexByScanNumberMap = new HashMap<String, Integer>();
			rtByScanNumberMap = new HashMap<String, Double>();
			BufferedReader br = null;
			int charge = -1;
			int scan1Num = -1;
			int scan2Num = -1;
			int numScan = -1;
			Double rt = null;
			try {
				br = new BufferedReader(new FileReader(ms2File));
				String line;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("S")) {
						numScan++;
						final String[] split = line.split("\\s");
						String scan1 = split[1];
						scan1Num = Integer.valueOf(scan1);
						String scan2 = split[2];
						scan2Num = Integer.valueOf(scan2);
						charge = -1;
					} else if (line.startsWith("Z")) {
						charge = Integer.valueOf(line.split("\\s")[1]);
						String key = scan1Num + "." + scan2Num + "." + charge;
						spectrumIndexByScanNumberMap.put(key, numScan);
						if (rt != null) {
							rtByScanNumberMap.put(key, rt);
						}
					} else if (line.startsWith("I")) {
						final String[] split = line.split("\\s");
						if (split[1].equalsIgnoreCase("RetTime")) {
							try {
								rt = Double.valueOf(split[2]);

							} catch (NumberFormatException e) {
								log.warn(e);
								rt = null;
							}
						}
					}
				}
				log.info(spectrumIndexByScanNumberMap.size() + " spectra readed in MS2 file. "
						+ ms2File.getAbsolutePath());
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
		}

	}

	public Double getSpectrumRTByScan(String scanNumber) {
		process();
		if (rtByScanNumberMap.containsKey(scanNumber)) {
			return rtByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	public String getFileName() {
		return FilenameUtils.getBaseName(ms2File.getAbsolutePath());
	}
}
