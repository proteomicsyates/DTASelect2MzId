package edu.scripps.yates.dtaselect2mzid.util;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.proteored.miapeapi.cv.Accession;
import org.proteored.miapeapi.cv.ControlVocabularyManager;
import org.proteored.miapeapi.cv.ControlVocabularySet;
import org.proteored.miapeapi.cv.ControlVocabularyTerm;
import org.proteored.miapeapi.cv.LocalOboTestControlVocabularyManager;
import org.proteored.miapeapi.cv.NEWTOntology;
import org.proteored.miapeapi.cv.PSIMassSpectrometryOntology;
import org.proteored.miapeapi.cv.UNIMODOntology;
import org.proteored.miapeapi.cv.UnitOntology;
import org.proteored.miapeapi.cv.ms.ContactPositionMS;

import edu.scripps.yates.dbindex.model.AssignMass;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.Role;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.model.utils.ModelConstants;

public class DTASelect2MzIdUtil {
	private static final Logger log = Logger.getLogger(DTASelect2MzIdUtil.class);
	private static final String UO_ONTOLOGY_ID = "UO";
	private static final String NEWT_ONTOLOGY_ID = "NEWT";
	private static final String PSI_ONTOLOGY_ID = "PSI-MS";
	private static final String UNIMOD_ONTOLOGY_ID = "UNIMOD";
	private static ControlVocabularyManager ontology = null;
	private static Map<String, ControlVocabularyTerm> cvTermCache = new HashMap<String, ControlVocabularyTerm>();

	public static ControlVocabularyManager getOntologyManager() {
		if (ontology == null) {
			log.info("Loading Controlled Vocabularies from Ontologies...");
			ontology = new LocalOboTestControlVocabularyManager();
			log.info("Ontologies loaded...");
		}
		return ontology;
	}

	public static Param getCVParam(String accession, String name, String value, Cv cv, String unitAccession, Cv unitCv,
			String unitName) {
		final Param param = new Param();
		final CvParam cvParam = new CvParam();
		cvParam.setAccession(accession);
		cvParam.setCv(cv);
		cvParam.setName(name);
		cvParam.setUnitAccession(unitAccession);
		cvParam.setUnitCv(unitCv);
		cvParam.setUnitName(unitName);
		cvParam.setValue(value);
		param.setParam(cvParam);
		return param;
	}

	public static Param getCVParam(ControlVocabularyTerm cvTerm, String value) {
		return getCVParam(cvTerm.getTermAccession().toUpperCase(), cvTerm.getPreferredName(), value,
				getCvByName(cvTerm.getCVRef()), null, null, null);
	}

	private static Cv getCvByName(String cvRef) {
		String cvId = getCVIdFromCVRef(cvRef);
		final CvList cvList = getCVList();
		for (Cv cv : cvList.getCv()) {
			if (cv.getId().equals(cvId)) {
				return cv;
			}
		}
		return null;
	}

	private static String getCVIdFromCVRef(String cvRef) {
		if ("MS".equals(cvRef)) {
			return getPSIMsCv().getId();
		}
		return cvRef;
	}

	public static CvList getCVList() {
		CvList ret = new CvList();
		ret.getCv().add(getPSIMsCv());
		ret.getCv().add(getNEWTCv());
		ret.getCv().add(getUoCv());
		ret.getCv().add(getUnimodCv());
		return ret;
	}

	public static Cv getUnimodCv() {
		Cv ret = new Cv();
		ret.setId(UNIMOD_ONTOLOGY_ID);
		ret.setFullName(UNIMODOntology.getFullName());
		ret.setUri(UNIMODOntology.getAddress());
		// ret.setVersion(UNIMODOntology.getVersion());
		return ret;
	}

	public static Cv getUoCv() {
		Cv ret = new Cv();
		ret.setId(UO_ONTOLOGY_ID);
		ret.setFullName(UnitOntology.getFullName());
		ret.setUri(UnitOntology.getAddress());
		// ret.setVersion(UnitOntology.getVersion());
		return ret;
	}

	public static Cv getNEWTCv() {
		Cv ret = new Cv();
		ret.setId(NEWT_ONTOLOGY_ID);
		ret.setFullName(NEWTOntology.getFullName());
		ret.setUri(NEWTOntology.getAddress());
		// ret.setVersion(NEWTOntology.getVersion());
		return ret;
	}

	public static Cv getPSIMsCv() {
		Cv ret = new Cv();
		ret.setId(PSI_ONTOLOGY_ID);
		ret.setFullName(PSIMassSpectrometryOntology.getFullName());
		ret.setUri(PSIMassSpectrometryOntology.getAddress());
		// ret.setVersion(PSIMassSpectrometryOntology.getVersion());
		return ret;
	}

	public static Param getCVParam(String accession, String name, String value, Cv cv) {
		return getCVParam(accession, name, value, cv, null, null, null);
	}

	public static Param getUserParam(String name, String value, Cv cv, String unitAccession, Cv unitCv,
			String unitName) {
		final Param param = new Param();
		final UserParam userParam = new UserParam();
		userParam.setName(name);
		userParam.setUnitAccession(unitAccession);
		userParam.setUnitCv(unitCv);
		userParam.setUnitName(unitName);
		userParam.setValue(value);
		param.setParam(userParam);
		return param;
	}

	public static Param getUserParam(String name, String value) {
		return getUserParam(name, value, null, null, null, null);
	}

	public static Role getResearcherRole() {
		Role ret = new Role();
		final String RESEARCHER_CV = "MS:1001271";
		final ControlVocabularyTerm cvTermByAcc = getCvTermByAcc(RESEARCHER_CV,
				ContactPositionMS.getInstance(getOntologyManager()));
		ret.setCvParam(getCVParam(cvTermByAcc, null).getCvParam());
		return ret;
	}

	public static Role getSoftwareVendorRole(String value) {
		Role ret = new Role();
		final String SOFT_VENDOR_CV = "MS:1001267";
		final ControlVocabularyTerm cvTermByAcc = getCvTermByAcc(SOFT_VENDOR_CV,
				ContactPositionMS.getInstance(getOntologyManager()));
		ret.setCvParam(getCVParam(cvTermByAcc, value).getCvParam());
		return ret;
	}

	public static ControlVocabularyTerm getCvTermByAcc(String acc, ControlVocabularySet cvSet) {
		if (cvTermCache.containsKey(acc)) {
			return cvTermCache.get(acc);
		}
		final ControlVocabularyTerm cvTermByAccession = getOntologyManager().getCVTermByAccession(new Accession(acc),
				cvSet);
		cvTermCache.put(acc, cvTermByAccession);
		if (cvTermByAccession == null) {
			log.warn(acc + " accession is not found in ontology.");
		}
		return cvTermByAccession;
	}

	public static String createMzIdentMLStartTag(String id, String version) {
		StringBuffer sb = new StringBuffer();
		// tag opening plus id attribute
		sb.append("<MzIdentML id=\"").append(id).append("\"");
		// further attributes
		sb.append(" version=\"" + version + "\"");
		if (version.equals("1.2.0")) {

			sb.append(" xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.2\"");
			sb.append(" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"");
			sb.append(
					" xsi:schemaLocation=\"http://psidev.info/psi/pi/mzIdentML/1.2 https://github.com/HUPO-PSI/mzIdentML/raw/master/schema/mzIdentML1.2.0-candidate.xsd\"");
			DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");
			sb.append(" creationDate=\"").append(dfm.format(new Date())).append("\"");
			// finally close the tag
			sb.append(" >");
		} else {
			sb.append(" xmlns=\"").append(ModelConstants.MZIDML_NS).append("\"");
			sb.append(" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"");
			sb.append(" xsi:schemaLocation=\"").append(ModelConstants.MZIDML_NS).append(" ")
					.append(ModelConstants.MZIDML_SCHEMA).append("\"");
			DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");
			sb.append(" creationDate=\"").append(dfm.format(new Date())).append("\"");
			// finally close the tag
			sb.append(" >");
		}
		return sb.toString();
	}

	public static Set<PSM> getPsms(Collection<Protein> dtaSelectProteins) {
		Set<PSM> psms = new HashSet<PSM>();
		for (Protein protein : dtaSelectProteins) {
			psms.addAll(protein.getPSMs());
		}
		return psms;
	}

	public static double getMonoMass(String sequence, ModificationParams modificationParams) {
		AssignMass.getInstance(true);
		double calculateMonoMass = 0.0;
		for (int index = 0; index < sequence.length(); index++) {
			final char aa = sequence.charAt(index);
			calculateMonoMass += AssignMass.getMass(aa);
			for (SearchModification searchModification : modificationParams.getSearchModification()) {
				if (searchModification.isFixedMod() && searchModification.getResidues().contains(aa)) {
					calculateMonoMass += searchModification.getMassDelta();
				}
			}
		}
		return calculateMonoMass;
	}

}
