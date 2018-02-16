
/**
 * <p>
 * Title:
 * </p>
 *
 * <p>
 * Description:
 * </p>
 *
 * <p>
 * Copyright: Copyright (c) 2008
 * </p>
 *
 * <p>
 * Company: Integrated Proteomics Solutions
 * </p>
 *
 * @author Tao Xu
 * @version $Id
 */

package edu.scripps.yates.dtaselect2mzid;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;

import org.apache.log4j.Logger;

import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Database;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.EnzymeInfo;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.EnzymeInfo.Residues;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.HeavySearchMode;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Isotopes;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.CTerm;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.CTerm.DiffMods;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.NTerm;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.NTerm.DiffMods.DiffMod;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.StaticMods;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.StaticMods.StaticMod;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.NumPeakLimits;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.PeptideLengthLimits;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.PrecursorMassLimits;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.SearchMode;
import edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Tolerance;

public class SearchXmlFile {
	private final static Logger log = Logger.getLogger(SearchXmlFile.class);
	private static JAXBContext jaxbContext;

	private Integer proteinDbId;
	private Long proteinUserId;
	private String databaseName = "";
	private String databasePath = "";
	private String primaryScoreType = "1";
	private String secondaryScoreType = "2";
	private String minMatch = "0";
	private String precursorIsotope = "mono";
	private String fragmentIsotope = "mono";
	private String numIsotopicPeaks = "3";
	private String preprocess = "1"; // preprocess mode
	private String highPrecursorTolerance = "3000";
	private String lowPrecursorTolerance = "3000";
	private String fragmentTolerance = "600";
	private String precursorTolerance = "50";
	private String maxPrecursorMass = "6000";
	private String minPrecursorMass = "600";
	private final String maxPrecursorCharge = "1000";
	private String minPrecursorCharge = "0";
	private String minimumPeptideLength = "6";
	private String maxNumPeaks = "5000";
	private String minNumPeaks = "25";
	// private String minNumPeaks = "8";
	private String maxNumDiffmod = "2";
	private String isDatabaseIndexed = "false";
	private String peakRankThreshold = "200";
	private String candidatePeptideThreshold = "500";
	private String numOutput = "5";
	private String locusType = "1";
	private String displayMod = "0";
	private String msaMode = "0";
	private String isDeCharged = "0";
	private String atomicEnrichment = "0";
	private String heavyAtomicEnrichment = "100";
	private String activationMethod = "CID";
	private final String paramFileName = "search.xml";
	private String chargeDisambiguation = "0";
	private String enzymeName = "Trypsin";
	private String enzymeSpecificity = "2";
	private String maxInternalMisCleavage = "-1"; // -1 for unlimited
	private String enzymeResidues = "";
	private String cutPosition = "true";
	private String[] ntermdiff;
	private String[] ctermdiff;
	private String[] diffmods;
	private String[] staticmods;
	private String includeHeavySearch = "no";
	private String mediumChecked = "no";
	private String isN15Search = "no";
	private String[] heavymods;
	private String[] mediummods;

	private String nTermStaticMod = "0";
	private String cTermStaticMod = "0";
	private String resolution = "low";
	private String[] resNames;
	private String pulseLightLabeling = null;
	private String pulseHeavyLabeling = null;

	private int scanNumber;
	private boolean exists;

	// Document object representing this search.xml file
	// private Document doc = null;

	public String getParamFileName() {
		return paramFileName;
	}

	public String getDatabasePath() {
		return databasePath;
	}

	public void setDatabasePath(String databasePath) {
		this.databasePath = databasePath;
	}

	public boolean isModificationSearch() {
		boolean withntermdiff = false;

		if (ntermdiff != null && ntermdiff.length > 0) {
			for (int i = 0; i < ntermdiff.length; i++) {
				if (ntermdiff[i] != null && !"".equals(ntermdiff[i].trim())) {
					String mass = ntermdiff[i].split("\\s+")[0];

					if ((!"".equals(mass)) && Double.parseDouble(mass) != 0)
						withntermdiff = true;
				}
			}
		}
		boolean withctermdiff = false;
		if (ctermdiff != null && ctermdiff.length > 0) {
			for (int i = 0; i < ctermdiff.length; i++) {
				if (ctermdiff[i] != null && !"".equals(ctermdiff[i].trim())) {
					String mass = ctermdiff[i].split("\\s+")[0];

					if ((!"".equals(mass)) && Double.parseDouble(mass) != 0)
						withctermdiff = true;
				}

			}
		}

		boolean withdiff = Integer.parseInt(maxNumDiffmod) > 0 ? true : false;

		return withntermdiff || withctermdiff || withdiff;
	}

	private void readExistingFile(String file) throws IOException, JAXBException {

		// System.out.println("Path: " + path + "\tFile: " + file );

		File paraFile = new File(file);

		if (!paraFile.exists()) {
			scanNumber = 4000;
			staticmods = new String[1];
			staticmods[0] = "57.02146 C";
			exists = false;
			return;
		}
		exists = true;
		if (jaxbContext == null) {
			try {
				jaxbContext = JAXBContext.newInstance(Parameters.class);
			} catch (JAXBException e) {
				e.printStackTrace();

			}
		}
		final Unmarshaller createUnmarshaller = jaxbContext.createUnmarshaller();
		final Parameters sp = (Parameters) createUnmarshaller.unmarshal(paraFile);

		SearchMode searchMode = sp.getSearchMode();
		if (searchMode != null) {
			primaryScoreType = searchMode.getPrimaryScoreType() + "";
			secondaryScoreType = searchMode.getSecondaryScoreType() + "";
			minMatch = searchMode.getMinMatch() + "";
			preprocess = searchMode.getPreprocess() + "";
			msaMode = searchMode.getMultistageActivationMode() + "";
			isDeCharged = String.valueOf(searchMode.getIsDecharged());
			peakRankThreshold = searchMode.getPeakRankThreshold() + "";
			candidatePeptideThreshold = searchMode.getCandidatePeptideThreshold() + "";
			numOutput = searchMode.getNumOutput() + "";
			locusType = searchMode.getLocusType() + "";
			chargeDisambiguation = String.valueOf(searchMode.getChargeDisambiguation());
			if (searchMode.getFragmentationMethod() != null) {
				activationMethod = searchMode.getFragmentationMethod().trim();
			}
		}
		Isotopes isotopes = sp.getIsotopes();
		if (isotopes != null) {
			precursorIsotope = isotopes.getPrecursor();
			fragmentIsotope = isotopes.getFragment();
			// numIsotopicPeaks = sp.getNumIsotopicPeaks() + "";
			String nips = Integer.valueOf(isotopes.getNumPeaks() + "") + "";
			numIsotopicPeaks = "0".equals(nips) ? "3" : nips;
			resolution = isotopes.getNumPeaks() == 0 ? "low" : "high";
		}
		Tolerance tolerance = sp.getTolerance();
		if (tolerance != null) {
			highPrecursorTolerance = tolerance.getPrecursorHigh() + "";
			lowPrecursorTolerance = tolerance.getPrecursorLow() + "";
			fragmentTolerance = tolerance.getFragment() + "";
			// System.out.println("fragmentTolerance: " + fragmentTolerance);
			precursorTolerance = tolerance.getPrecursor() + "";
		}
		PrecursorMassLimits precursorMassLimits = sp.getPrecursorMassLimits();
		if (precursorMassLimits != null) {
			maxPrecursorMass = precursorMassLimits.getMaximum() + "";
			minPrecursorMass = precursorMassLimits.getMinimum() + "";
		}
		PeptideLengthLimits peptideLengthLimits = sp.getPeptideLengthLimits();
		if (peptideLengthLimits != null) {
			minimumPeptideLength = peptideLengthLimits.getMinimum() + "";
		}
		NumPeakLimits numPeakLimits = sp.getNumPeakLimits();
		if (numPeakLimits != null) {
			maxNumPeaks = numPeakLimits.getMaximum() + "";
			minNumPeaks = numPeakLimits.getMinimum() + "";
		}
		maxNumDiffmod = sp.getMaxNumDiffmod() + "";
		Database database = sp.getDatabase();
		if (database != null) {
			databaseName = database.getDatabaseName();
			isDatabaseIndexed = database.getIsIndexed();
		}

		Modifications modifications = sp.getModifications();
		if (modifications != null) {
			displayMod = modifications.getDisplayMod() + "";
			NTerm nTerm = modifications.getNTerm();
			if (nTerm != null) {
				edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.NTerm.StaticMod nTermStaticModObj = nTerm
						.getStaticMod();
				if (nTermStaticModObj != null) {
					nTermStaticMod = nTermStaticModObj.getMassShift() + "";

				}
				edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.NTerm.DiffMods nTermDiffModsObj = nTerm
						.getDiffMods();
				if (nTermDiffModsObj != null) {
					List<DiffMod> nTermDiffModList = nTermDiffModsObj.getDiffMod();
					if (nTermDiffModList != null) {
						ntermdiff = new String[nTermDiffModList.size()];
						int counter = 0;
						for (DiffMod nTermDiffMod : nTermDiffModList) {
							ntermdiff[counter++] = nTermDiffMod.getMassShift() + "";
						}
					}
				}
			}
			CTerm cTerm = modifications.getCTerm();
			if (cTerm != null) {
				edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.CTerm.StaticMod cTermStaticModObj = cTerm
						.getStaticMod();
				if (cTermStaticModObj != null) {
					cTermStaticMod = cTermStaticModObj.getMassShift() + "";
				}
				DiffMods cTermDiffModObj = cTerm.getDiffMods();
				if (cTermDiffModObj != null) {
					List<edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.CTerm.DiffMods.DiffMod> cTermDiffModList = cTermDiffModObj
							.getDiffMod();
					if (cTermDiffModList != null) {
						ctermdiff = new String[cTermDiffModList.size()];
						int counter = 0;
						for (edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.CTerm.DiffMods.DiffMod cTermDiffMod : cTermDiffModList) {
							ctermdiff[counter++] = cTermDiffMod.getMassShift() + "";
						}

					}
				}
			}
			edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.DiffMods diffModsObj = modifications
					.getDiffMods();
			if (diffModsObj != null) {
				List<edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.DiffMods.DiffMod> diffModList = diffModsObj
						.getDiffMod();
				if (diffModList != null) {
					diffmods = new String[diffModList.size()];
					int counter = 0;
					for (edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.DiffMods.DiffMod diffMod : diffModList) {
						StringBuffer sb = new StringBuffer(20);
						sb.append(diffMod.getMassShift() + " ");
						edu.scripps.yates.dtaselect2mzid.search.jaxb.Parameters.Modifications.DiffMods.DiffMod.Residues diffModResiduesObj = diffMod
								.getResidues();
						if (diffModResiduesObj != null) {
							final List<String> diffModResiduesList = diffModResiduesObj.getResidue();
							for (String diffModResidue : diffModResiduesList) {
								sb.append(diffModResidue);
							}
						}
						diffmods[counter++] = sb.toString();
					}
				}
			}
			// get static mods
			StaticMods staticModObj = modifications.getStaticMods();
			List<StaticMod> staticModList = staticModObj.getStaticMod();
			if (staticModList != null) {
				staticmods = new String[staticModList.size()];
				int counter = 0;
				for (StaticMod staticMod : staticModList) {
					char residue = 0;
					if (staticMod.getResidue() != null) {
						residue = staticMod.getResidue().charAt(0);
					}
					String massShift = staticMod.getMassShift() + "";
					staticmods[counter++] = massShift + " " + residue;
				}
			}
		}
		// paramFileName = "search.xml";

		EnzymeInfo enzymeInfo = sp.getEnzymeInfo();
		if (enzymeInfo != null) {
			enzymeSpecificity = enzymeInfo.getSpecificity() + "";
			maxInternalMisCleavage = enzymeInfo.getMaxNumInternalMisCleavage() + "";
			// System.out.println("sp:" + sp);
			// System.out.println("protease: " + sp.getProtease());
			// System.out.println("residues: " +
			// sp.getProtease().getResidues());
			Residues residues = enzymeInfo.getResidues();
			if (residues != null) {
				for (String residue : residues.getResidue()) {
					enzymeResidues += residue;
				}
			}
			// cutPosition = sp.getProtease().getType() ? "true" : "false";
			enzymeName = enzymeInfo.getName();
		}

		FileInputStream paramInput = new FileInputStream(paraFile);

		// root.getChild("database");
		// root.getChild("search_mode");
		// root.getChild("isotopes");
		// root.getChild("tolerance");
		// root.getChild("precursor_mass_limits");
		// atomicEnrichment = sp.getAtomicEnrichment() + "";
		// activationMethod = sp.isEtdSearch()? "ETD" : "CID";

		HeavySearchMode hsm = sp.getHeavySearchMode();
		if (hsm != null) {
			includeHeavySearch = "yes";
			isN15Search = hsm.getIsN15Search();

			if ("yes".equals(isN15Search)) {
				String haenrich = hsm.getHeavyAtomicEnrichment() + "";
				if (haenrich != null) {
					setHeavyAtomicEnrichment(haenrich);
				}
			}
			List<String> heavyModList = hsm.getHeavyMods().getHeavyMod();
			int numheavymods = heavyModList.size();
			// System.out.println("Numheavymods: " + numheavymods);
			if (numheavymods > 0) {
				ArrayList<String> allhmods = new ArrayList<String>(20);
				for (String hmod : heavyModList) {
					// System.out.println("heavymod: " + e.getTextTrim() +
					// "\tindex: " + counter);
					// need to check the

					// logger.debug("hmod: " + hmod);
					if (hmod != null && !hmod.equals("")) {
						allhmods.add(hmod);
					}

				}

				int numhmods = allhmods.size() > 0 ? allhmods.size() : 1;
				heavymods = new String[numhmods];
				int counter = 0;
				if (numhmods == 1) {
					if (allhmods.size() > 0) {

						heavymods[0] = allhmods.get(0);
					} else {
						heavymods[0] = heavyModList.get(0);

					}
				} else {
					for (String s : allhmods) {
						heavymods[counter++] = s;
					}

				}
			}
			// for mediumList
			// List<Element> mediumModList =
			// hsm.getChild("medium_mods").getChildren();
			// int numMediummods = mediumModList.size();
			List<String> mediumModList = new ArrayList<String>();
			int numMediummods = 0;
			if (hsm.getMediumMods() != null) {
				mediumModList = hsm.getMediumMods().getMediumMod();
				numMediummods = mediumModList.size();
			}

			// System.out.println("numMediummods: " + numMediummods);
			if (numMediummods > 0) {
				ArrayList<String> allmmods = new ArrayList<String>(20);
				for (String mmod : mediumModList) {
					// String mmod = e.getTextTrim();
					if (mmod != null && !mmod.equals("")) {
						allmmods.add(mmod);
					}
				}
				int nummods = allmmods.size() > 0 ? allmmods.size() : 1;
				mediummods = new String[nummods];
				int counter = 0;
				if (nummods == 1) {
					if (allmmods.size() > 0)
						mediummods[0] = allmmods.get(0);
					else
						mediummods[0] = mediumModList.get(0);
				} else {
					for (String s : allmmods)
						mediummods[counter++] = s;
				}
			}
		}

		// "peptide_length_limits");
		// private String includeHeavySearch = "no";
		// private String isN15Search = "no";
		// private String [] heavymods;
		paramInput.close();

	}

	public SearchXmlFile(String filename) throws IOException, JAXBException {
		readExistingFile(filename);
	}

	public void setDiffmods(String[] d) {
		diffmods = d;
	}

	public String[] getDiffmods() {
		return diffmods;
	}

	public int getNumDiffmodType() {
		if (diffmods != null)
			return diffmods.length;
		else
			return 0;
	}

	public void setStaticmods(String[] d) {
		staticmods = d;
	}

	public void setHeavymods(String[] d) {
		heavymods = d;
	}

	public String[] getHeavymods() {
		return heavymods;
	}

	public void setResNames(String[] d) {
		resNames = d;
	}

	public String[] getResNames() {
		return resNames;
	}

	public void setMediummods(String[] d) {
		mediummods = d;
	}

	public String[] getMediummods() {
		return mediummods;
	}

	public String[] getStaticmods() {
		return staticmods;
	}

	public void setNtermdiff(String[] d) {
		ntermdiff = d;
	}

	public String[] getNtermdiff() {
		return ntermdiff;
	}

	public void setCtermdiff(String[] d) {
		ctermdiff = d;
	}

	public String[] getCtermdiff() {
		return ctermdiff;
	}

	public void setNTermStaticMod(String s) {
		nTermStaticMod = s;
	}

	public String getNTermStaticMod() {
		return nTermStaticMod;
	}

	public void setCTermStaticMod(String s) {
		cTermStaticMod = s;
	}

	public String getCTermStaticMod() {
		return cTermStaticMod;
	}

	public void setNumOutput(String o) {
		numOutput = o;
	}

	public String getNumOutput() {
		return numOutput;
	}

	public void setAtomicEnrichment(String c) {
		atomicEnrichment = c;
	}

	public String getAtomicEnrichment() {
		return atomicEnrichment;
	}

	public void setChargeDisambiguation(String c) {
		chargeDisambiguation = c;
	}

	public String getChargeDisambiguation() {
		return chargeDisambiguation;
	}

	public void setCutPosition(String e) {
		cutPosition = e;
	}

	public String getCutPosition() {
		return cutPosition;
	}

	public void setResolution(String r) {
		resolution = r;
	}

	public String getResolution() {
		return resolution;
	}

	public void setEnzymeName(String e) {
		enzymeName = e;
	}

	public String getEnzymeName() {
		return enzymeName;
	}

	public void setEnzymeSpecificity(String e) {
		enzymeSpecificity = e;
	}

	public String getEnzymeSpecificity() {
		return enzymeSpecificity;
	}

	public void setMaxInternalMisCleavage(String msc) {
		maxInternalMisCleavage = msc;
	}

	public String getMaxInternalMisCleavage() {
		return maxInternalMisCleavage;
	}

	public void setEnzymeResidues(String r) {
		enzymeResidues = r;

	}

	public String getEnzymeResidues() {
		return enzymeResidues;
	}

	public void setMaxNumPeaks(String numSpectra) {
		maxNumPeaks = numSpectra;
	}

	public String getMaxNumPeaks() {
		return maxNumPeaks;
	}

	public String setMinNumPeaks() {
		return minNumPeaks;
	}

	public void setMinNumPeaks(String numSpectra) {
		minNumPeaks = numSpectra;
	}

	public String getMaxNumDiffmod() {
		return maxNumDiffmod;
	}

	public void setMaxNumDiffmod(String max) {
		if (max == null || "".equals(max.trim())) {
			max = "0";
		}
		if (max.indexOf(".") > 0) {
			max = max.split("\\.")[0];
		}
		maxNumDiffmod = max;
	}

	public String getLocusType() {
		return locusType;
	}

	public Integer getProteinDbId() {
		return proteinDbId;
	}

	public void setProteinDbId(Integer id) {
		proteinDbId = id;
	}

	public Long getProteinUserId() {
		return proteinUserId;
	}

	public void setProteinUserId(Long proteinUserId) {
		this.proteinUserId = proteinUserId;
	}

	public String getIsDeCharged() {
		return isDeCharged;
	}

	public void setIsDeCharged(String c) {
		isDeCharged = c;
	}

	public String getActivationMethod() {
		return activationMethod;
	}

	public void setActivationMethod(String c) {
		activationMethod = c;
	}

	public void setPeakRankThreshold(String p) {
		peakRankThreshold = p;
	}

	public String getPeakRankThreshold() {
		return peakRankThreshold;
	}

	public void setCandidatePeptideThreshold(String s) {
		candidatePeptideThreshold = s;
	}

	public String getCandidatePeptideThreshold() {
		return candidatePeptideThreshold;
	}

	public String getIsDatabaseIndexed() {
		return isDatabaseIndexed;
	}

	public void setIsDatabaseIndexed(String i) {
		isDatabaseIndexed = i;
	}

	public String getPrimaryScoreType() {
		return primaryScoreType;
	}

	public void setPrimaryScoreType(String type) {
		primaryScoreType = type;
	}

	public String getSecondaryScoreType() {
		return secondaryScoreType;
	}

	public void setSecondaryScoreType(String type) {
		secondaryScoreType = type;
	}

	public String getMinMatch() {
		return minMatch;
	}

	public void setMinMatch(String m) {
		minMatch = m;
	}

	public void setNumIsotopicPeaks(String s) {
		numIsotopicPeaks = s;
	}

	public String getNumIsotopicPeaks() {
		return numIsotopicPeaks;
	}

	public String getPrecursorIsotope() {
		return precursorIsotope;
	}

	public void setPrecursorIsotope(String p) {
		precursorIsotope = p;
	}

	public String getFragmentIsotope() {
		return fragmentIsotope;
	}

	public void setFragmentIsotope(String f) {
		fragmentIsotope = f;
	}

	public void setPreprocess(String c) {
		preprocess = c;
	}

	public String getPreprocess() {
		return preprocess;
	}

	public void setHighPrecursorTolerance(String s) {
		highPrecursorTolerance = s;
		lowPrecursorTolerance = s;
	}

	public String getHighPrecursorTolerance() {
		return highPrecursorTolerance;
	}

	public void setLowPrecursorTolerance(String s) {
		lowPrecursorTolerance = s;
		highPrecursorTolerance = s;
	}

	public String getLowPrecursorTolerance() {
		return lowPrecursorTolerance;
	}

	public String getPrecursorTolerance() {
		return precursorTolerance;
	}

	public void setPrecursorTolerance(String pt) {
		precursorTolerance = pt;
	}

	public String getFragmentTolerance() {
		return fragmentTolerance;
	}

	public void setFragmentTolerance(String ft) {
		fragmentTolerance = ft;
	}

	public String getMinimumPeptideLength() {
		return minimumPeptideLength;
	}

	public void setMinimumPeptideLength(String m) {
		minimumPeptideLength = m;
	}

	public void setMaxPrecursorMass(String s) {
		maxPrecursorMass = s;
	}

	public String getMaxPrecursorMass() {
		return maxPrecursorMass;
	}

	public void setMinPrecursorMass(String s) {
		minPrecursorMass = s;
	}

	public String getMinPrecursorMass() {
		return minPrecursorMass;
	}

	public void setMinPrecursorCharge(String s) {
		minPrecursorCharge = s;
	}

	public String getMinPrecursorCharge() {
		return minPrecursorCharge;
	}

	public void setMsaMode(String s) {
		msaMode = s;
	}

	public String getMsaMode() {
		return msaMode;
	}

	public void setDisplayMod(String s) {
		displayMod = s;
	}

	public String getDisplayMod() {
		return displayMod;
	}

	public String getMinNumPeaks() {
		return minNumPeaks;
	}

	public String getDatabaseName() {
		return databaseName;
	}

	public String getDatabaseNameWithoutPath() {

		int splitIndex = databaseName.lastIndexOf(File.separator);
		if (splitIndex <= 0)
			return databaseName;

		return databaseName.substring(databaseName.lastIndexOf(File.separator));
	}

	public void setDatabaseName(String dbname) {
		databaseName = dbname;
	}

	public void setIncludeHeavySearch(String b) {
		includeHeavySearch = b;
	}

	public String getIncludeHeavySearch() {
		return includeHeavySearch;
	}

	public void setMediumChecked(String mc) {
		mediumChecked = mc;
	}

	public String getMediumChecked() {
		return mediumChecked;
	}

	public void setIsN15Search(String b) {
		isN15Search = b;
	}

	public String getIsN15Search() {
		return isN15Search;
	}

	public void setHeavyAtomicEnrichment(String b) {
		heavyAtomicEnrichment = b;
	}

	public String getHeavyAtomicEnrichment() {
		return heavyAtomicEnrichment;
	}

	public String output() {
		StringBuffer sb = new StringBuffer(10000);
		sb.append("<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n");
		sb.append("<!--  Parameters for ProLuCID Seach  -->\n");
		sb.append("<parameters>\n");
		sb.append("\t<database>\n");
		sb.append("\t\t<database_name>" + getDatabaseName() + "</database_name>\n");
		sb.append("\t\t<is_indexed>" + getIsDatabaseIndexed() + "<is_indexed>");
		sb.append("\t</database>\n");

		sb.append("</parameters>\n");
		return sb.toString();
	}

	public void outputSequestParamsWithLocalPath(String fullpath) throws FileNotFoundException {
		PrintStream ps = new PrintStream(fullpath);
		// logger.debug("In outputSequestParamsWithLocalPath, databasePath is "
		// + getDatabasePath() + "\t and databaseNameWithoutPath is " +
		// getDatabaseNameWithoutPath() + "\tdatabaseName is "
		// +getDatabaseName() );
		ps.println("[SEQUEST]");
		ps.println("database_name = " + getDatabasePath() + File.separator + getDatabaseNameWithoutPath());
		ps.close();
	}

	public void outputFakeSequestParams(String fullpath) throws FileNotFoundException {
		PrintStream ps = new PrintStream(fullpath);
		// logger.debug("In outputFakeSequestParams, databasePath is " +
		// getDatabasePath() + "\t and databaseName is " + getDatabaseName() + "
		// and database is outputed wihtout adding path");
		ps.println("[SEQUEST]");
		ps.println("database_name = " + getDatabaseName());
		ps.close();
	}

	public static boolean isHeavySearchFolder(String searchfolder) throws IOException {
		String searchxmlfile = searchfolder + "/search.xml";
		String sequestparamfile = searchfolder + "/sequest.params";
		File searchparam = new File(searchxmlfile);
		File sequestparam = new File(sequestparamfile);
		boolean isheavysearch = false;
		if (searchparam.exists()) {
			BufferedReader br = null;
			try {
				br = new BufferedReader(new FileReader(searchparam));
				String line = "";
				while ((line = br.readLine()) != null) {
					if (line.contains("heavy_search_mode")) {
						isheavysearch = true;
						break;
					}
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			} finally {
				if (br != null)
					br.close();
			}

		} else if (sequestparam.exists()) {

			BufferedReader br = null;
			try {
				br = new BufferedReader(new FileReader(sequestparam));
				String line = "";
				while ((line = br.readLine()) != null) {
					if (line.contains("heavy_search")) {
						isheavysearch = true;
						break;
					}
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			} finally {
				if (br != null)
					br.close();
			}
		}

		return isheavysearch;
	}

	public static void main(String args[]) throws Exception {
		// SearchXmlFile sp = new SearchXmlFile();
		// System.out.println(sp.output());
		// sp.outputXml(System.out);
		/*
		 * String hfolder =
		 * "/data/2/rpark/ip2_data/mesingh/HAPTestRatBrainJun2010/UP_HAP1_2010_06_28_14_1434/search/projects2010_08_05_134984/";
		 * String lfolder =
		 * "/data/2/rpark/ip2_data/taoxu/abrf2010/step12_2010_01_15_18_560/search/projects2010_10_17_116432";
		 * System.out.println("Is heavy search in " + hfolder + ": " +
		 * isHeavySearchFolder(hfolder)); System.out.println(
		 * "Is heavy search in " + lfolder + ": " +
		 * isHeavySearchFolder(lfolder));
		 */

		SearchXmlFile sp = new SearchXmlFile(
				"/data/2/rpark/ip2_data/taoxu/abrf2010/step12_2010_01_15_18_560/search/projects2010_10_17_116432/");

		System.out.println(sp.getFragmentTolerance());
	}

	public String getPulseLightLabeling() {
		return pulseLightLabeling;
	}

	public void setPulseLightLabeling(String pulseLightLabeling) {
		this.pulseLightLabeling = pulseLightLabeling;
	}

	public String getPulseHeavyLabeling() {
		return pulseHeavyLabeling;
	}

	public void setPulseHeavyLabeling(String pulseHeavyLabeling) {
		this.pulseHeavyLabeling = pulseHeavyLabeling;
	}

	public int getScanNumber() {
		return scanNumber;
	}

	public void setScanNumber(int scanNumber) {
		this.scanNumber = scanNumber;
	}

	public boolean exists() {
		return exists;
	}

}
