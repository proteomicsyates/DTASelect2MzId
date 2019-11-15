package edu.scripps.yates.dtaselect2mzid.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class MS2Reader {
	private final static Logger log = Logger.getLogger(MS2Reader.class);
	private final File ms2File;
	private TObjectIntHashMap<String> spectrumIndexByScanNumberMap;
	private TObjectDoubleHashMap<String> rtByScanNumberMap;
	private TObjectDoubleHashMap<String> precursorIntensityByScanNumberMap;
	private TObjectDoubleHashMap<String> precursorMassByScanNumberMap;

	private boolean grabPrecursorIntensitiesAndMass;

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
			spectrumIndexByScanNumberMap = new TObjectIntHashMap<String>();
			precursorIntensityByScanNumberMap = new TObjectDoubleHashMap<String>();
			precursorMassByScanNumberMap = new TObjectDoubleHashMap<String>();
			rtByScanNumberMap = new TObjectDoubleHashMap<String>();
			BufferedReader br = null;
			int charge = -1;
			int scan1Num = -1;
			int scan2Num = -1;
			int numScan = -1;
			Double rt = null;
			Double intensity = null;
			Double mass = null;
			try {
				br = new BufferedReader(new FileReader(ms2File));
				String line;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("S")) {
						numScan++;
						final String[] split = line.split("\\s");
						final String scan1 = split[1];
						scan1Num = Integer.valueOf(scan1);
						final String scan2 = split[2];
						scan2Num = Integer.valueOf(scan2);
						charge = -1;
						if (grabPrecursorIntensitiesAndMass && split.length > 3) {
							mass = Double.valueOf(split[3]);
						}
					} else if (line.startsWith("Z")) {
						charge = Integer.valueOf(line.split("\\s")[1]);
						final String key = scan1Num + "." + scan2Num + "." + charge;
						spectrumIndexByScanNumberMap.put(key, numScan);
						if (rt != null) {
							rtByScanNumberMap.put(key, rt);
						}
						if (intensity != null) {
							precursorIntensityByScanNumberMap.put(key, intensity);
						}
						if (mass != null) {
							precursorMassByScanNumberMap.put(key, mass);
						}
					} else if (line.startsWith("I")) {
						final String[] split = line.split("\\s");
						if (split[1].equalsIgnoreCase("RetTime") || split[1].equalsIgnoreCase("RTime")) {
							try {
								rt = Double.valueOf(split[2]);

							} catch (final NumberFormatException e) {
								log.warn(e);
								rt = null;
							}
						} else if (grabPrecursorIntensitiesAndMass) {
							if (split[1].equalsIgnoreCase("PrecursorInt")) {
								try {
									intensity = Double.valueOf(split[2]);

								} catch (final NumberFormatException e) {
									log.warn(e);
									intensity = null;
								}
							}
						}
					}
				}
				log.info(spectrumIndexByScanNumberMap.size() + " spectra read in MS2 file. "
						+ ms2File.getAbsolutePath());
			} catch (final FileNotFoundException e) {
				e.printStackTrace();
			} catch (final IOException e) {
				e.printStackTrace();
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (final IOException e) {
						e.printStackTrace();
					}
				}
			}
		}

	}

	/**
	 *
	 * @param scanNumber in the form of "scan.scan.charge"
	 * @return
	 */
	public Double getSpectrumRTByScan(String scanNumber) {
		process();
		if (rtByScanNumberMap.containsKey(scanNumber)) {
			return rtByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	/**
	 *
	 * @param scanNumber in the form of "scan.scan.charge"
	 * @return
	 */
	public Double getPrecursorIntensityByScan(String scanNumber) {
		process();
		if (precursorIntensityByScanNumberMap.containsKey(scanNumber)) {
			return precursorIntensityByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	/**
	 *
	 * @param scanNumber in the form of "scan.scan.charge"
	 * @return
	 */
	public Double getPrecursorMassByScan(String scanNumber) {
		process();
		if (precursorMassByScanNumberMap.containsKey(scanNumber)) {
			return precursorMassByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	public String getFileName() {
		return FilenameUtils.getBaseName(ms2File.getAbsolutePath());
	}

	public TObjectIntHashMap<String> getSpectrumIndexByScanNumberMap() {
		if (spectrumIndexByScanNumberMap == null) {
			process();
		}
		return spectrumIndexByScanNumberMap;
	}

	public TObjectDoubleHashMap<String> getRtByScanNumberMap() {
		if (rtByScanNumberMap == null) {
			process();
		}
		return rtByScanNumberMap;
	}

	public boolean isGrabPrecursorIntensities() {
		return grabPrecursorIntensitiesAndMass;
	}

	public void setGrabPrecursorIntensities(boolean grabPrecursorIntensities) {
		this.grabPrecursorIntensitiesAndMass = grabPrecursorIntensities;
	}

}
