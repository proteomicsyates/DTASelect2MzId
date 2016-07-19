package edu.scripps.yates.dtaselect2mzid;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.JAXBException;

import org.proteored.miapeapi.cv.Accession;
import org.proteored.miapeapi.cv.ControlVocabularyManager;
import org.proteored.miapeapi.cv.ControlVocabularyTerm;
import org.proteored.miapeapi.cv.msi.CleavageName;

import com.compomics.dbtoolkit.io.EnzymeLoader;

import edu.scripps.yates.dtaselect2mzid.util.DTASelect2MzIdUtil;
import edu.scripps.yates.dtaselect2mzid.util.LabeledSearchType;
import edu.scripps.yates.dtaselect2mzid.util.PeptideModificationUtil;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.Enzymes;
import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
import uk.ac.ebi.jmzidml.model.mzidml.ParamList;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;

public class SearchParametersManager {
	/** modification params by search type **/
	private final Map<LabeledSearchType, ModificationParams> modificationParamsByLST = new HashMap<LabeledSearchType, ModificationParams>();
	private final Map<LabeledSearchType, SearchXmlFile> searchXmlFileByLST = new HashMap<LabeledSearchType, SearchXmlFile>();
	private static File inputFolder;
	private static SearchParametersManager instance;
	private final ControlVocabularyManager cvManager = DTASelect2MzIdUtil.getOntologyManager();
	private final Map<LabeledSearchType, Boolean> someDiffModByLst = new HashMap<LabeledSearchType, Boolean>();
	private final Map<LabeledSearchType, Boolean> someFixedModByLst = new HashMap<LabeledSearchType, Boolean>();

	public static SearchParametersManager getInstance() throws IOException {
		if (inputFolder == null) {
			throw new IOException(
					"no input folder is configured. Call to the static method setInputFolder before getting any instance");
		}
		if (instance == null) {
			instance = new SearchParametersManager();
		}
		return instance;
	}

	private SearchParametersManager() {

	}

	public static void setInputFolder(File inputFolder) {
		SearchParametersManager.inputFolder = inputFolder;
	}

	public SearchXmlFile getSearchParameters(LabeledSearchType lst) {
		if (!searchXmlFileByLST.containsKey(lst)) {
			try {
				final String searchXMLFileName = inputFolder.getAbsolutePath() + File.separator + lst.getKey()
						+ "search.xml";
				SearchXmlFile searchXmlFile = new SearchXmlFile(searchXMLFileName);
				if (searchXmlFile.exists()) {
					searchXmlFileByLST.put(lst, searchXmlFile);
				} else {
					searchXmlFileByLST.put(lst, null);
				}
			} catch (IOException e) {
				e.printStackTrace();
			} catch (JAXBException e) {
				e.printStackTrace();
			}
		}
		return searchXmlFileByLST.get(lst);
	}

	public ModificationParams getModificationParams(LabeledSearchType lst) {
		final SearchXmlFile searchParameters = getSearchParameters(lst);
		if (searchParameters == null) {
			return null;
		}
		if (!modificationParamsByLST.containsKey(lst)) {
			boolean containsValidMod = false;
			ModificationParams modificationParams = new ModificationParams();
			boolean someDiffMod = false;
			boolean someFixedMod = false;
			// N term fix
			if (searchParameters.getNTermStaticMod() != null) {

				final Float massDiff = Float.valueOf(searchParameters.getNTermStaticMod());
				if (Float.compare(massDiff, 0.0f) != 0) {
					containsValidMod = true;
					someFixedMod = true;
					modificationParams.getSearchModification()
							.add(getSearchModification(".", massDiff, true, true, false));
				}
			}
			// N term
			if (searchParameters.getNtermdiff() != null) {
				Set<String> set = new HashSet<String>();
				for (String nTermDiff : searchParameters.getNtermdiff()) {
					if (!set.contains(nTermDiff)) {
						set.add(nTermDiff);
						final Float massDiff = Float.valueOf(nTermDiff);
						if (Float.compare(massDiff, 0.0f) != 0) {
							containsValidMod = true;
							someDiffMod = true;
							modificationParams.getSearchModification()
									.add(getSearchModification(".", massDiff, false, true, false));
						}
					}
				}
			}
			// C term fix
			if (searchParameters.getCTermStaticMod() != null) {
				final Float massDiff = Float.valueOf(searchParameters.getCTermStaticMod());
				if (Float.compare(massDiff, 0.0f) != 0) {
					containsValidMod = true;
					someFixedMod = true;
					modificationParams.getSearchModification()
							.add(getSearchModification(".", massDiff, true, false, true));
				}
			}
			// C term
			if (searchParameters.getCtermdiff() != null) {
				containsValidMod = true;
				Set<String> set = new HashSet<String>();
				for (String cTermDiff : searchParameters.getCtermdiff()) {
					if (!set.contains(cTermDiff)) {
						set.add(cTermDiff);
						final Float massDiff = Float.valueOf(cTermDiff);
						if (Float.compare(massDiff, 0.0f) != 0) {
							containsValidMod = true;
							someDiffMod = true;
							modificationParams.getSearchModification()
									.add(getSearchModification(".", massDiff, false, false, true));
						}
					}
				}
			}
			// modifications fix
			if (searchParameters.getStaticmods() != null) {
				Set<String> set = new HashSet<String>();
				for (String diffAndResidue : searchParameters.getStaticmods()) {
					if (!set.contains(diffAndResidue)) {
						containsValidMod = true;
						someFixedMod = true;
						set.add(diffAndResidue);
						Float massDiff = Float.valueOf(diffAndResidue.split(" ")[0]);
						String residues = diffAndResidue.split(" ")[1];
						modificationParams.getSearchModification()
								.add(getSearchModification(residues, massDiff, true, false, false));
					}
				}
			}
			// modifications
			if (searchParameters.getDiffmods() != null) {
				for (String diffAndResidue : searchParameters.getDiffmods()) {
					containsValidMod = true;
					someDiffMod = true;
					Float massDiff = Float.valueOf(diffAndResidue.split(" ")[0]);
					String residues = diffAndResidue.split(" ")[1];
					modificationParams.getSearchModification()
							.add(getSearchModification(residues, massDiff, false, false, false));
				}
			}

			someDiffModByLst.put(lst, someDiffMod);
			someFixedModByLst.put(lst, someFixedMod);

			if (containsValidMod) {
				return modificationParams;
			}
		}
		return modificationParamsByLST.get(lst);
	}

	ParamList getAdditionalSearchParams(LabeledSearchType lst, int numFileIDs) {
		ParamList ret = new ParamList();
		final SearchXmlFile searchParameters = getSearchParameters(lst);
		if (searchParameters == null) {
			return null;
		}

		if (searchParameters.getActivationMethod() != null && !"".equals(searchParameters.getActivationMethod())) {
			if ("CID".equalsIgnoreCase(searchParameters.getActivationMethod())) {
				ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1000133", "collision-induced dissociation", null,
						DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
			} else if ("ETD".equalsIgnoreCase(searchParameters.getActivationMethod())) {
				ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1000598", "electron transfer dissociation", null,
						DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
			} else {
				ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1000044", "Activation Method",
						searchParameters.getActivationMethod(), DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
			}
		}
		if (searchParameters.getPrimaryScoreType() != null && !"".equals(searchParameters.getPrimaryScoreType())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("primary_score_type", searchParameters.getPrimaryScoreType()).getUserParam());
		}
		if (searchParameters.getSecondaryScoreType() != null && !"".equals(searchParameters.getSecondaryScoreType())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("secondary_score_type", searchParameters.getSecondaryScoreType()).getUserParam());
		}
		if (searchParameters.getLocusType() != null && !"".equals(searchParameters.getLocusType())) {
			ret.getUserParam()
					.add(DTASelect2MzIdUtil.getUserParam("locus_type", searchParameters.getLocusType()).getUserParam());
		}
		if (searchParameters.getChargeDisambiguation() != null
				&& !"".equals(searchParameters.getChargeDisambiguation())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("charge_disambiguation", searchParameters.getChargeDisambiguation()).getUserParam());
		}
		if (searchParameters.getAtomicEnrichment() != null && !"".equals(searchParameters.getAtomicEnrichment())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("atomic_enrichement", searchParameters.getAtomicEnrichment()).getUserParam());
		}
		if (searchParameters.getMinMatch() != null && !"".equals(searchParameters.getMinMatch())) {
			ret.getUserParam()
					.add(DTASelect2MzIdUtil.getUserParam("min_match", searchParameters.getMinMatch()).getUserParam());
		}
		if (searchParameters.getPeakRankThreshold() != null && !"".equals(searchParameters.getPeakRankThreshold())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("peak_rank_threshold", searchParameters.getPeakRankThreshold()).getUserParam());
		}
		if (searchParameters.getCandidatePeptideThreshold() != null
				&& !"".equals(searchParameters.getCandidatePeptideThreshold())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("candidate_peptide_threshold", searchParameters.getCandidatePeptideThreshold())
					.getUserParam());
		}
		if (searchParameters.getCandidatePeptideThreshold() != null
				&& !"".equals(searchParameters.getCandidatePeptideThreshold())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("candidate_peptide_threshold", searchParameters.getCandidatePeptideThreshold())
					.getUserParam());
		}
		if (searchParameters.getNumOutput() != null && !"".equals(searchParameters.getNumOutput())) {
			ret.getUserParam()
					.add(DTASelect2MzIdUtil.getUserParam("num_output", searchParameters.getNumOutput()).getUserParam());
		}
		if (searchParameters.getCandidatePeptideThreshold() != null
				&& !"".equals(searchParameters.getCandidatePeptideThreshold())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("candidate_peptide_threshold", searchParameters.getCandidatePeptideThreshold())
					.getUserParam());
		}
		if (searchParameters.getIsDeCharged() != null && !"".equals(searchParameters.getIsDeCharged())) {
			ret.getUserParam().add(
					DTASelect2MzIdUtil.getUserParam("is_decharged", searchParameters.getIsDeCharged()).getUserParam());
		}
		if (searchParameters.getPreprocess() != null && !"".equals(searchParameters.getPreprocess())) {
			ret.getUserParam().add(
					DTASelect2MzIdUtil.getUserParam("preprocess", searchParameters.getPreprocess()).getUserParam());
		}
		if (searchParameters.getMinPrecursorMass() != null && !"".equals(searchParameters.getMinPrecursorMass())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("minimum_precursor_mass", searchParameters.getMinPrecursorMass()).getUserParam());
		}
		if (searchParameters.getMaxPrecursorMass() != null && !"".equals(searchParameters.getMaxPrecursorMass())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("maximum_precursor_mass", searchParameters.getMaxPrecursorMass()).getUserParam());
		}
		if (searchParameters.getMinPrecursorCharge() != null && !"".equals(searchParameters.getMinPrecursorCharge())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("minimum_precursor_charge", searchParameters.getMinPrecursorCharge()).getUserParam());
		}
		if (searchParameters.getMaxPrecursorMass() != null && !"".equals(searchParameters.getMaxPrecursorMass())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("maximum_precursor_charge", searchParameters.getMaxPrecursorMass()).getUserParam());
		}
		if (searchParameters.getMinimumPeptideLength() != null
				&& !"".equals(searchParameters.getMinimumPeptideLength())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("minimum_peptide_length", searchParameters.getMinimumPeptideLength()).getUserParam());
		}
		if (searchParameters.getMinNumPeaks() != null && !"".equals(searchParameters.getMinNumPeaks())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("minimum_num_peaks", searchParameters.getMinNumPeaks()).getUserParam());
		}
		if (searchParameters.getMaxNumPeaks() != null && !"".equals(searchParameters.getMaxNumPeaks())) {
			ret.getUserParam().add(DTASelect2MzIdUtil
					.getUserParam("maximum_num_peaks", searchParameters.getMaxNumPeaks()).getUserParam());
		}
		if (searchParameters.getMaxNumDiffmod() != null && !"".equals(searchParameters.getMaxNumDiffmod())) {
			ret.getUserParam()
					.add(DTASelect2MzIdUtil
							.getUserParam("maximum_num_differential_modifications", searchParameters.getMaxNumDiffmod())
							.getUserParam());
		}
		final String fragmentIsotope = searchParameters.getFragmentIsotope();
		if ("mono".equals(fragmentIsotope)) {
			ret.getCvParam()
					.add(DTASelect2MzIdUtil
							.getCVParam("MS:1001256", "fragment mass type mono", null, DTASelect2MzIdUtil.getPSIMsCv())
							.getCvParam());
		} else if ("avg".equals(fragmentIsotope)) {
			ret.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1001255", "fragment mass type average", null, DTASelect2MzIdUtil.getPSIMsCv())
					.getCvParam());
		}
		final String precursorIsotope = searchParameters.getPrecursorIsotope();
		if ("mono".equals(precursorIsotope)) {
			ret.getCvParam()
					.add(DTASelect2MzIdUtil
							.getCVParam("MS:1001211", "parent mass type mono", null, DTASelect2MzIdUtil.getPSIMsCv())
							.getCvParam());
		} else if ("avg".equals(precursorIsotope)) {
			ret.getCvParam()
					.add(DTASelect2MzIdUtil
							.getCVParam("MS:1001212", "parent mass type average", null, DTASelect2MzIdUtil.getPSIMsCv())
							.getCvParam());
		}

		// sample prefractionation:
		if (numFileIDs > 1) {
			ret.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1002493", "sample pre-fractionation", "MudPIT", DTASelect2MzIdUtil.getPSIMsCv())
					.getCvParam());
		}
		return ret;
	}

	public Tolerance getParentTolerance(LabeledSearchType lst) {
		Tolerance ret = new Tolerance();
		String plusValue = null;
		String minusValue = null;
		final SearchXmlFile searchParameters = getSearchParameters(lst);
		if (searchParameters == null) {
			return null;
		}
		if (searchParameters.getPrecursorIsotope().equals("mono")) {
			plusValue = searchParameters.getPrecursorTolerance();
			minusValue = searchParameters.getPrecursorTolerance();
		} else {
			plusValue = searchParameters.getHighPrecursorTolerance();
			minusValue = searchParameters.getLowPrecursorTolerance();
		}
		CvParam fragmentCVPlus = DTASelect2MzIdUtil.getCVParam("MS:1001412", "search tolerance plus value", plusValue,
				DTASelect2MzIdUtil.getPSIMsCv(), "UO:0000169", DTASelect2MzIdUtil.getUoCv(), "parts per million")
				.getCvParam();
		ret.getCvParam().add(fragmentCVPlus);
		CvParam fragmentCVMinus = DTASelect2MzIdUtil
				.getCVParam("MS:1001413", "search tolerance minus value", minusValue, DTASelect2MzIdUtil.getPSIMsCv(),
						"UO:0000169", DTASelect2MzIdUtil.getUoCv(), "parts per million")
				.getCvParam();
		ret.getCvParam().add(fragmentCVMinus);
		return ret;
	}

	public Tolerance getFragmentTolerance(LabeledSearchType lst) {
		Tolerance ret = new Tolerance();
		final SearchXmlFile searchParameters = getSearchParameters(lst);
		if (searchParameters == null) {
			return null;
		}
		CvParam fragmentCVPlus = DTASelect2MzIdUtil.getCVParam("MS:1001412", "search tolerance plus value",
				searchParameters.getFragmentTolerance(), DTASelect2MzIdUtil.getPSIMsCv(), "UO:0000169",
				DTASelect2MzIdUtil.getUoCv(), "parts per million").getCvParam();
		ret.getCvParam().add(fragmentCVPlus);
		CvParam fragmentCVMinus = DTASelect2MzIdUtil.getCVParam("MS:1001413", "search tolerance minus value",
				searchParameters.getFragmentTolerance(), DTASelect2MzIdUtil.getPSIMsCv(), "UO:0000169",
				DTASelect2MzIdUtil.getUoCv(), "parts per million").getCvParam();
		ret.getCvParam().add(fragmentCVMinus);
		return ret;
	}

	public Enzymes getEnzymes(LabeledSearchType lst) {
		Enzymes ret = new Enzymes();
		final Enzyme enzyme = getEnzyme(lst);
		if (enzyme == null) {
			return null;
		}
		ret.getEnzyme().add(enzyme);
		return ret;
	}

	private Enzyme getEnzyme(LabeledSearchType lst) {
		Enzyme ret = new Enzyme();
		final SearchXmlFile searchParameters = getSearchParameters(lst);
		if (searchParameters == null) {
			return null;
		}
		ret.setMissedCleavages(Integer.valueOf(searchParameters.getMaxInternalMisCleavage()));
		String enzymeName = searchParameters.getEnzymeName();
		if (enzymeName.equalsIgnoreCase("LysC")) {
			enzymeName = "Lys-C";
		}
		final Accession acc = cvManager.getControlVocabularyId(enzymeName, CleavageName.getInstance(cvManager));
		if (acc != null) {
			final ControlVocabularyTerm cvTerm = cvManager.getCVTermByAccession(acc,
					CleavageName.getInstance(cvManager));
			ParamList paramList = new ParamList();
			paramList.getCvParam().add(DTASelect2MzIdUtil.getCVParam(cvTerm, null).getCvParam());
			ret.setEnzymeName(paramList);
			ret.setId("ENZYME_" + cvTerm.getPreferredName() + "_" + lst.getKey());

		} else {
			ParamList paramList = new ParamList();
			paramList.getUserParam()
					.add(DTASelect2MzIdUtil.getUserParam(searchParameters.getEnzymeName(), null).getUserParam());
			ret.setEnzymeName(paramList);
			ret.setId("ENZYME_" + searchParameters.getEnzymeName() + "_" + lst.getKey());
		}
		if (searchParameters.getEnzymeResidues() != null && !"".equals(searchParameters.getEnzymeResidues())) {
			ret.setSiteRegexp(searchParameters.getEnzymeResidues());
		} else {
			try {
				final com.compomics.util.protein.Enzyme enzyme = EnzymeLoader.loadEnzyme(enzymeName, null);
				if (enzyme.getRestrict() != null) {
					StringBuilder sb = new StringBuilder();
					for (char aa : enzyme.getCleavage()) {
						if (!"".equals(sb.toString())) {
							sb.append("|");
						}
						sb.append(aa);
					}
					ret.setSiteRegexp(sb.toString());
				}
			} catch (IOException e) {
			}
		}
		return ret;

	}

	private SearchModification getSearchModification(String residues, float massDelta, boolean fixed, boolean nterm,
			boolean cterm) {
		SearchModification ret = new SearchModification();
		ret.setMassDelta(massDelta);
		if (residues != null) {
			for (int i = 0; i < residues.length(); i++) {
				String aa = String.valueOf(residues.charAt(i));
				ret.getResidues().add(aa);
			}
		}
		ret.setFixedMod(fixed);
		PeptideModificationUtil peptideModUtil = new PeptideModificationUtil(massDelta, residues);
		if (peptideModUtil.getAccession() != null) {
			ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam(peptideModUtil.getAccession(), peptideModUtil.getName(),
					null, DTASelect2MzIdUtil.getUnimodCv()).getCvParam());
		} else {
			ret.getCvParam()
					.add(DTASelect2MzIdUtil
							.getCVParam("MS:1001460", "unknown modification", null, DTASelect2MzIdUtil.getUnimodCv())
							.getCvParam());
		}
		if (nterm) {
			ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1001189", "modification specificity peptide N-term",
					null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		}
		if (cterm) {
			ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1001190", "modification specificity peptide C-term",
					null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		}
		return ret;
	}
}
