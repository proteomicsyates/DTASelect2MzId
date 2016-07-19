package edu.scripps.yates.dtaselect2mzid;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.proteored.miapeapi.cv.ControlVocabularyManager;
import org.proteored.miapeapi.cv.ControlVocabularyTerm;
import org.proteored.miapeapi.cv.ms.SoftwareName;
import org.proteored.miapeapi.cv.msi.Score;
import org.proteored.miapeapi.cv.msi.SearchType;

import edu.scripps.yates.dbindex.model.AssignMass;
import edu.scripps.yates.dtaselect.ProteinDTASelectParser;
import edu.scripps.yates.dtaselect.ProteinImplFromDTASelect;
import edu.scripps.yates.dtaselect2mzid.util.DTASelect2MzIdUtil;
import edu.scripps.yates.dtaselect2mzid.util.LabeledSearchType;
import edu.scripps.yates.dtaselect2mzid.util.MS2Reader;
import edu.scripps.yates.dtaselect2mzid.util.MzIdentMLVersion;
import edu.scripps.yates.dtaselect2mzid.util.PeptideModificationUtil;
import edu.scripps.yates.dtaselect2mzid.util.PeptideResultsDetailsCVSet;
import edu.scripps.yates.dtaselect2mzid.util.ReferenceToSpectra;
import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.utilities.checksum.MD5Checksum;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.masses.MassesUtil;
import edu.scripps.yates.utilities.model.factories.ProteinEx;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.util.Pair;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSampleCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.ContactRole;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.FileFormat;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentationTable;
import uk.ac.ebi.jmzidml.model.mzidml.InputSpectra;
import uk.ac.ebi.jmzidml.model.mzidml.InputSpectrumIdentifications;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.MassTable;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
import uk.ac.ebi.jmzidml.model.mzidml.Organization;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.ParamList;
import uk.ac.ebi.jmzidml.model.mzidml.ParentOrganization;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinAmbiguityGroup;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetection;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.Residue;
import uk.ac.ebi.jmzidml.model.mzidml.Role;
import uk.ac.ebi.jmzidml.model.mzidml.Sample;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabase;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabaseRef;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SourceFile;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIDFormat;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentification;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItemRef;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;

public class DTASelect2MzId {
	private static final Logger log = Logger.getLogger(DTASelect2MzIdUtil.class);
	private static final String SEQUEST_CV = "MS:1001208";
	// TODO include it when supported by skyline
	private static final String PROLUCID_CV = "MS:1002596";
	private static final String PROTEIN_DESCRIPTION_CV = "MS:1001088";
	private static final String DTASELECT_FORMAT_CV = "MS:1001464";
	private static final String DTASELECT_TOOL_CV = "MS:1002598";
	private static final String MS2_FORMAT_CV = "MS:1001466";
	private static final String MSMSSearchType = "MS:1001083";

	private static final String YATES_LAB_CONTACT_ID = "YATES_LAB";
	private static final String SCRIPPS_ORG_ID = "TSRI";
	private static final String DELTACN_SEQUEST = "MS:1001156";
	private static final String DELTACN_PROLUCID = "MS:1002535";
	private static final String XCORR_SEQUEST = "MS:1001155";
	private static final String XCORR_PROLUCID = "MS:1002534";
	// private static final String DTASELECT_CV = "MS:1002598";
	private static final String SCAN_START_TIME_CV = "MS:1000016";
	private static final String DISTINCT_PEP_SEQUENCSE = "MS:1001097";
	private static final String PROTEIN_COVERAGE_CV = "MS:1001093";
	private static final String EMPAI_VALUE_CV = "MS:1001905";
	private static final String NSAF_CV = "MS:1002547";
	private static final String LEADING_PROTEIN = "MS:1002401";
	private static final String NON_LEADING_PROTEIN = "MS:1002402";
	private static final String NO_FIXED_MODIFICATIONS_SEARCHED = "MS:1002453";
	private static final String NO_VARIABLE_MODIFICATIONS_SEARCHED = "MS:1002454";
	// TODO REMOVE
	private static final String FAKE_SCORE_FOR_SKYLINE = "MS:1002467";
	private static final String SEQUENCE_SAME_SET_PROTEIN = "MS:1001594";
	private static final String SEQUENCE_SUBSET_SET_PROTEIN = "MS:1001596";
	private static final String FINAL_PSM_LIST_CV = "MS:1002439";
	private static final String MzXML_FORMAT_CV = "MS:1000566";
	private static final String COUNT_IDENTIFIED_PROTEINS_CV = "MS:1002404";
	private static final String PG_PASS_THRESHOLD_CV = "MS:1002415";
	private static Options options;
	private final List<File> dtaSelectFiles = new ArrayList<File>();
	private final File output;
	private ProteinDTASelectParser dtaSelectParser;

	private SearchDatabase searchDatabase;

	private final Map<String, DBSequence> dbSequences = new HashMap<String, DBSequence>();
	private final Map<String, Peptide> peptides = new HashMap<String, Peptide>();
	private final Map<String, PeptideEvidence> peptideEvidencesByKey = new HashMap<String, PeptideEvidence>();
	private final Map<String, MS2Reader> ms2ReaderByFileName = new HashMap<String, MS2Reader>();
	private final Map<String, SpectrumIdentificationList> spectrumIdentificationListMap = new HashMap<String, SpectrumIdentificationList>();
	private final Map<String, SpectrumIdentificationResult> sirMap = new HashMap<String, SpectrumIdentificationResult>();
	private final Map<String, SpectrumIdentificationItem> siiMap = new HashMap<String, SpectrumIdentificationItem>();
	private Sample sample;
	private final ControlVocabularyManager cvManager = DTASelect2MzIdUtil.getOntologyManager();
	private final Map<String, ProteinAmbiguityGroup> pag = new HashMap<String, ProteinAmbiguityGroup>();
	private final Map<String, SpectraData> spectraDataBySpectraFileName = new HashMap<String, SpectraData>();
	private String decoyPatternString;
	private final Map<LabeledSearchType, SpectrumIdentificationProtocol> sipByLST = new HashMap<LabeledSearchType, SpectrumIdentificationProtocol>();
	private ProteinDetectionProtocol pdp;
	private List<ProteinGroup> groups;
	private boolean ignoreSpectra = false;
	private Pattern decoyPattern;
	private final DecimalFormat myFormatter = new DecimalFormat("#.##");
	private final DecimalFormat myFormatter4digits = new DecimalFormat("#.####");
	private final Map<String, SpectrumIdentification> siMap = new HashMap<String, SpectrumIdentification>();
	private final Map<String, PeptideHypothesis> phByPeptideEvidence = new HashMap<String, PeptideHypothesis>();
	private final Map<String, ProteinDetectionHypothesis> pdhMapByID = new HashMap<String, ProteinDetectionHypothesis>();
	private final Map<String, LabeledSearchType> labeledSearchTypeByFileName = new HashMap<String, LabeledSearchType>();
	private final Map<ProteinEx, Set<Protein>> proteinSetByProteinEx = new HashMap<ProteinEx, Set<Protein>>();
	private ProteinDetectionList pdl;
	private ReferenceToSpectra referenceToSpectra = ReferenceToSpectra.MS2;

	// by default, version 1_1
	private final MzIdentMLVersion mzIdentMLVersion;
	private boolean skylineCompatible;

	private static final Set<String> spectraFileNamesNotFound = new HashSet<String>();;

	public DTASelect2MzId(Collection<File> dtaSelectFiles, File outputMzIdentMLFile, MzIdentMLVersion version) {
		if (dtaSelectFiles == null) {
			throw new IllegalArgumentException("DTASelect files are null");
		}
		for (File file : dtaSelectFiles) {
			addDTASelect(file);
		}
		output = outputMzIdentMLFile;
		mzIdentMLVersion = version;
	}

	public DTASelect2MzId(File dtaSelectFile, File outputMzIdentMLFile, MzIdentMLVersion version) {
		if (dtaSelectFile == null) {
			throw new IllegalArgumentException("DTASelect file is not found");
		}

		addDTASelect(dtaSelectFile);
		output = outputMzIdentMLFile;
		mzIdentMLVersion = version;
	}

	private void addDTASelect(File file) {
		if (!file.exists()) {
			throw new IllegalArgumentException("File " + FilenameUtils.getName(file.getAbsolutePath())
					+ " is not found at '" + file.getParentFile().getAbsolutePath() + "'");
		}
		if (!dtaSelectFiles.contains(file)) {
			dtaSelectFiles.add(file);
		}
	}

	public void setReferenceToSpectra(ReferenceToSpectra referenceToSpectra) {
		this.referenceToSpectra = referenceToSpectra;
	}

	/**
	 * mzIdentML<br>
	 * cvList<br>
	 * AnalysisSoftwareList<br>
	 * Provider<br>
	 * AuditCollection<br>
	 * AnalysisSampleCollection<br>
	 * SequenceCollection<br>
	 * AnalysisCollection<br>
	 * AnalysisProtocolCollection<br>
	 * DataCollection<br>
	 * Inputs<br>
	 * AnalysisData<br>
	 * SpectrumIdentificationList<br>
	 * ProteinDetectionList<br>
	 * /AnalysisData<br>
	 * /DataCollection<br>
	 * BibliographicReference<br>
	 * /mzIdentML<br>
	 *
	 * @throws IOException
	 */
	public void convert() throws IOException {
		SearchParametersManager.setInputFolder(getInputFolder());

		dtaSelectParser = new ProteinDTASelectParser(dtaSelectFiles);
		dtaSelectParser.setDecoyPattern(decoyPatternString);

		if (!ignoreSpectra && referenceToSpectra == ReferenceToSpectra.MS2) {
			final String path = getDTASelectFolder(dtaSelectFiles);
			final Set<String> spectraFileFullPaths = dtaSelectParser.getSpectraFileFullPaths();
			for (String spectraFileFullPath : spectraFileFullPaths) {
				// MS2Reader
				try {
					File ms2File = new File(spectraFileFullPath);
					if (!ms2File.exists()) {
						ms2File = new File(path + File.separator + spectraFileFullPath + ".ms2");
					}
					MS2Reader ms2Reader = new MS2Reader(ms2File);
					log.info("Setting up ms2 reader for file: '" + ms2File.getAbsolutePath());
					String spectraFileName = FilenameUtils.getBaseName(ms2File.getAbsolutePath());
					if (ms2ReaderByFileName.containsKey(spectraFileName)) {
						log.info(ms2ReaderByFileName.get(spectraFileName).getFileName());
					}

					if (getLabeledSearchTypeByFileName(spectraFileName) != LabeledSearchType.LIGHT) {
						spectraFileName = spectraFileName.substring(1, spectraFileName.length());
					}
					ms2ReaderByFileName.put(spectraFileName, ms2Reader);
				} catch (IllegalArgumentException e) {
					log.error(e.getMessage());
				}
			}
		}
		// Note: writing of '\n' characters is optional and only for readability
		// of the produced XML document
		// Also note: since the XML is produced in individual parts, the overall
		// formatting of the document
		// is not as nice as it would be when marshalling the whole structure at
		// once.
		FileWriter writer = null;
		try {
			writer = new FileWriter(output);
			MzIdentMLMarshaller m = new MzIdentMLMarshaller();
			// XML header
			log.info("Creating XML header...");
			writer.write(m.createXmlHeader() + "\n");
			// mzIdentML start tag
			writer.write(createMzIdentMLStartTag(m, mzIdentMLVersion, "DTASelect2MzId") + "\n");

			log.info("Creating CvList element...");
			CvList cvList = DTASelect2MzIdUtil.getCVList();
			m.marshal(cvList, writer);
			writer.write("\n");

			//
			log.info("Creating AnalysisSoftwareList element...");
			AnalysisSoftwareList analysisSoftwareList = getAnalysisSoftwareList();
			m.marshal(analysisSoftwareList, writer);
			writer.write("\n");

			log.info("Creating Provider element...");
			Provider provider = getDefaultProvider();
			m.marshal(provider, writer);
			writer.write("\n");

			//
			//
			log.info("Creating AuditCollection element...");
			AuditCollection auditCollection = getAuditCollection();
			m.marshal(auditCollection, writer);
			writer.write("\n");
			//
			log.info("Creating AnalysisSampleCollection element...");
			AnalysisSampleCollection analysisSampleCollection = getAnalysisSampleCollection(true);
			if (analysisSampleCollection != null) {
				m.marshal(analysisSampleCollection, writer);
				writer.write("\n");
			}
			//
			log.info("Creating SequenceCollection element...");
			SequenceCollection sequenceCollection = getSequenceCollection();
			log.info("Sequence collection with " + sequenceCollection.getDBSequence().size() + " DBSequences, "
					+ sequenceCollection.getPeptide().size() + " Peptides and "
					+ sequenceCollection.getPeptideEvidence().size() + " PeptideEvidences");
			m.marshal(sequenceCollection, writer);
			writer.write("\n");

			log.info("Creating SpectrumIdentificationResults...");
			// get fisrt the spectra
			for (PSM psm : dtaSelectParser.getPSMsByPSMID().values()) {
				// System.out.println(psm.getMSRun().getRunId() + "\t" +
				// psm.getPSMIdentifier());
				SpectrumIdentificationResult specIdentRes = getSpectrumIdentificationResult(psm);
			}
			log.info(dtaSelectParser.getPSMsByPSMID().size() + " SpectrumIdentificationResult elements created");

			//
			log.info("Creating AnalysisCollection element...");
			AnalysisCollection analysisCollection = getAnalysisCollection();
			m.marshal(analysisCollection, writer);
			writer.write("\n");
			//
			log.info("Creating AnalysisProcoloCollection element...");
			AnalysisProtocolCollection analysisProtocolCollection = getAnalysisProcotolCollection();
			m.marshal(analysisProtocolCollection, writer);
			writer.write("\n");
			//
			log.info("Creating DataCollection element...");
			writer.write(m.createDataCollectionStartTag() + "\n");
			//
			log.info("Creating Inputs element...");
			Inputs inputs = getInputs();
			m.marshal(inputs, writer);
			writer.write("\n");
			//
			log.info("Creating AnalysisData element...");
			writer.write(m.createAnalysisDataStartTag() + "\n");
			//
			// one spectrumIdentificationList by each fileName
			// (pre-fractionation) and each labelledSearchTag if
			// available
			for (String fileID : dtaSelectParser.getSpectraFileNames()) {
				for (LabeledSearchType lst : LabeledSearchType.values()) {
					final SpectrumIdentificationList sil = getSpectrumIdentificationList(fileID, lst);
					if (sil != null) {

						m.marshal(sil, writer);
						writer.write("\n");

					} else {
						// log.info("Search '" + fileID + "' - '" + lst.name() +
						// "' not found.");
					}
				}
			}
			//
			log.info("Creating ProteinDetectionList element...");
			final ProteinDetectionList proteinDetectionList = getProteinDetectionList();
			// writer.write(m.createProteinDetectionListStartTag(proteinDetectionList.getId(),
			// null) + "\n");
			m.marshal(proteinDetectionList, writer);
			// // group with panalyzer
			// List<ProteinGroup> groups2 = getGroups();
			// for (ProteinGroup proteinGroup : groups2) {
			// ProteinAmbiguityGroup protAmbGroup =
			// getProteinAmbiguityGroup(proteinGroup);
			// m.marshal(protAmbGroup, writer);
			// writer.write("\n");
			// }
			log.info(groups.size() + " ProteinAmbiguityGroup elements created in ProteinDetectionList.");
			//
			// writer.write(m.createProteinDetectionListClosingTag() + "\n");
			//
			writer.write(m.createAnalysisDataClosingTag() + "\n");
			//
			writer.write(m.createDataCollectionClosingTag() + "\n");
			//
			// BibliographicReference ref =
			// unmarshaller.unmarshal(MzIdentMLElement.BibliographicReference.getXpath());
			// m.marshal(ref, writer);
			// writer.write("\n");

			writer.write(m.createMzIdentMLClosingTag());
			log.info("File created at: " + output.getAbsolutePath());

			final Set<String> errorMessages = PeptideModificationUtil.getErrorMessages();
			if (!errorMessages.isEmpty()) {
				System.err.println(
						"Some modifications were not identified properly and were translated as 'unknown' modifications:");
				for (String errorMessage : errorMessages) {
					System.err.println(errorMessage);
				}
			}
		} finally {
			if (writer != null)
				writer.close();
			if (mzIdentMLVersion == MzIdentMLVersion.VERSION_1_2) {
				fixLines();
			}
		}
	}

	private void fixLines() {
		try {
			Path path = Paths.get(output.toURI());
			Charset charset = StandardCharsets.UTF_8;

			String content = new String(Files.readAllBytes(path), charset);
			content = content.replaceAll("http://psidev.info/psi/pi/mzIdentML/1.1",
					"http://psidev.info/psi/pi/mzIdentML/1.2");
			Files.write(path, content.getBytes(charset));

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private String createMzIdentMLStartTag(MzIdentMLMarshaller marshaller, MzIdentMLVersion version, String id) {
		switch (version) {
		case VERSION_1_1:
			return marshaller.createMzIdentMLStartTag(id);
		case VERSION_1_2:
			return createMzIdentMLStartTagForVersion_1_2(id);
		default:
			throw new IllegalArgumentException(version.name() + " is not supported.");
		}
	}

	private String createMzIdentMLStartTagForVersion_1_2(String id) {
		StringBuffer sb = new StringBuffer();
		String MZIDML_NS = "http://psidev.info/psi/pi/mzIdentML/1.2";
		String MZIDML_VERSION = "1.2.0";
		String MZIDML_SCHEMA = "https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0-candidate.xsd";
		String MZIDSCHEMA_LOCATION = "http://psidev.info/psi/pi/mzIdentML/1.2 " + MZIDML_SCHEMA;
		// tag opening plus id attribute
		sb.append("<MzIdentML id=\"").append(id).append("\"");
		// further attributes
		sb.append(" version=\"").append(MZIDML_VERSION).append("\"");
		sb.append(" xmlns=\"").append(MZIDML_NS).append("\"");
		sb.append(" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"");
		sb.append(" xsi:schemaLocation=\"").append(MZIDSCHEMA_LOCATION).append("\"");
		DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");
		sb.append(" creationDate=\"").append(dfm.format(getCreationDate())).append("\"");
		// finally close the tag
		sb.append(" >");

		return sb.toString();
	}

	private Date getCreationDate() {
		if (dtaSelectFiles != null && !dtaSelectFiles.isEmpty()) {
			for (File file : dtaSelectFiles) {
				try {
					BasicFileAttributes attr = Files.readAttributes(file.toPath(), BasicFileAttributes.class);
					final Date date = new Date(attr.creationTime().toMillis());
					return date;
				} catch (IOException e) {
				}
			}
		}
		return new Date();

	}

	private String getDTASelectFolder(List<File> dtaSelectFiles2) {

		return dtaSelectFiles2.get(0).getParentFile().getAbsolutePath();
	}

	private List<ProteinGroup> getGroups() {
		if (groups == null) {
			log.info("Grouping  " + dtaSelectParser.getProteins().size()
					+ " proteins according to PAnalyzer algorithm...");
			groups = new ArrayList<ProteinGroup>();
			Set<GroupableProtein> groupableProteins = new HashSet<GroupableProtein>();
			for (Set<Protein> proteinSet : dtaSelectParser.getProteins().values()) {

				groupableProteins.addAll(proteinSet);
			}
			PAnalyzer panalyzer = new PAnalyzer(true);
			groups = panalyzer.run(groupableProteins);
			log.info("Proteins grouped in " + groups.size() + " groups.");
		}
		return groups;

	}

	private ProteinAmbiguityGroup getProteinAmbiguityGroup(ProteinGroup group) throws IOException {
		String pagKey = getProteinAmbiguityGroupKey(group);
		if (pag.containsKey(pagKey)) {
			return pag.get(pagKey);
		}
		ProteinAmbiguityGroup ret = new ProteinAmbiguityGroup();
		pag.put(pagKey, ret);
		ret.setId(pagKey);
		final Map<String, Set<Protein>> proteinsByAcc = getProteinsByAcc(group);
		final Map<String, ProteinEvidence> evidences = new HashMap<String, ProteinEvidence>();
		for (String acc : proteinsByAcc.keySet()) {
			evidences.put(acc, proteinsByAcc.get(acc).iterator().next().getEvidence());
		}
		for (String acc : proteinsByAcc.keySet()) {
			Set<Protein> proteinSet = proteinsByAcc.get(acc);
			ret.getProteinDetectionHypothesis().add(getProteinDetectionHypothesis(proteinSet));
		}

		ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam(PG_PASS_THRESHOLD_CV, "protein group passes threshold",
				"true", DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		return ret;
	}

	private Map<String, Set<Protein>> getProteinsByAcc(ProteinGroup group) {
		Map<String, Set<Protein>> map = new HashMap<String, Set<Protein>>();
		for (GroupableProtein groupableProtein : group) {
			final String accession = groupableProtein.getAccession();
			if (map.containsKey(accession)) {
				map.get(accession).add((Protein) groupableProtein);
			} else {
				Set<Protein> set = new HashSet<Protein>();
				set.add((Protein) groupableProtein);
				map.put(accession, set);
			}
		}
		return map;
	}

	private String getProteinAmbiguityGroupKey(ProteinGroup group) throws IOException {
		StringBuilder sb = new StringBuilder();
		List<String> pdhIds = new ArrayList<String>();
		final Map<String, Set<Protein>> proteinsByAcc = getProteinsByAcc(group);
		for (Set<Protein> proteinSet : proteinsByAcc.values()) {
			final String pdhId = getProteinDetectionHypothesis(proteinSet).getId();
			if (!pdhIds.contains(pdhId)) {
				pdhIds.add(pdhId);
			}
		}
		for (String pdhId : pdhIds) {
			if (!"".equals(sb.toString())) {
				sb.append("_");
			}
			sb.append(pdhId);
		}
		return "PAG_" + sb.toString();
	}

	/**
	 *
	 * @param dtaSelectProteins
	 *            can be more than one because can be the same protein from
	 *            different msRuns
	 * @return
	 * @throws IOException
	 */
	private ProteinDetectionHypothesis getProteinDetectionHypothesis(Collection<Protein> dtaSelectProteins)
			throws IOException {
		// we assume that the evidence of all proteins in the collection is the
		// same
		Protein protein1 = dtaSelectProteins.iterator().next();
		final DBSequence dbSequence = getDBSequence(protein1);
		final String pdhid = getPDHID(dbSequence);
		if (pdhMapByID.containsKey(pdhid)) {
			return pdhMapByID.get(pdhid);
		}
		ProteinDetectionHypothesis pdh = new ProteinDetectionHypothesis();
		pdhMapByID.put(pdhid, pdh);
		pdh.setId(pdhid);
		pdh.setDBSequence(dbSequence);
		pdh.setPassThreshold(true);
		final Cv psiMsCv = DTASelect2MzIdUtil.getPSIMsCv();
		pdh.getCvParam()
				.add(DTASelect2MzIdUtil.getCVParam("MS:1001592", "family member protein", null, psiMsCv).getCvParam());
		final ProteinEvidence evidence = protein1.getEvidence();
		switch (evidence) {
		case AMBIGUOUSGROUP:
			pdh.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1002216", "PAnalyzer:ambiguous group member", null, psiMsCv).getCvParam());
			pdh.getCvParam()
					.add(DTASelect2MzIdUtil.getCVParam(LEADING_PROTEIN, "leading protein", null, psiMsCv).getCvParam());
			break;
		case CONCLUSIVE:
			pdh.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1002213", "PAnalyzer:conclusive protein", null, psiMsCv).getCvParam());
			pdh.getCvParam()
					.add(DTASelect2MzIdUtil.getCVParam(LEADING_PROTEIN, "leading protein", null, psiMsCv).getCvParam());
			break;
		case FILTERED:
			break;
		case INDISTINGUISHABLE:
			pdh.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1002214", "PAnalyzer:indistinguishable protein", null, psiMsCv).getCvParam());
			pdh.getCvParam()
					.add(DTASelect2MzIdUtil.getCVParam(LEADING_PROTEIN, "leading protein", null, psiMsCv).getCvParam());
			pdh.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam(SEQUENCE_SAME_SET_PROTEIN, "sequence same-set protein", null, psiMsCv).getCvParam());

			break;
		case NONCONCLUSIVE:
			pdh.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1002215", "PAnalyzer:non-conclusive protein", null, psiMsCv).getCvParam());
			pdh.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam(NON_LEADING_PROTEIN, "non-leading protein", null, psiMsCv).getCvParam());
			pdh.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam(SEQUENCE_SUBSET_SET_PROTEIN, "sequence sub-set protein", null, psiMsCv).getCvParam());
			break;
		default:
			break;
		}

		// pdh.getCvParam().add(DTASelect2MzIdUtil
		// .getCVParam("MS:1001194", "quality estimation with decoy database",
		// null, psiMsCv).getCvParam());
		// if (getProteinThreshold() != null) {
		// final String thresholdValue =
		// getProteinThreshold().getCvParam().get(0).getValue();
		// pdh.getCvParam().add(DTASelect2MzIdUtil
		// .getCVParam("MS:1001447", "prot:FDR threshold", thresholdValue,
		// psiMsCv).getCvParam());
		// }

		if (protein1 instanceof ProteinImplFromDTASelect) {
			ProteinImplFromDTASelect proteinTmp = (ProteinImplFromDTASelect) protein1;
			// protein coverage
			pdh.getCvParam().add(DTASelect2MzIdUtil.getCVParam(PROTEIN_COVERAGE_CV, "sequence coverage",
					myFormatter.format(proteinTmp.getCoverage()), psiMsCv).getCvParam());
			// empai value
			if (proteinTmp.getEMPai() != null) {
				pdh.getCvParam().add(DTASelect2MzIdUtil
						.getCVParam(EMPAI_VALUE_CV, "emPAI value", myFormatter.format(proteinTmp.getEMPai()), psiMsCv)
						.getCvParam());
			}
			// NSAF
			pdh.getCvParam().add(DTASelect2MzIdUtil.getCVParam(NSAF_CV, "normalized spectral abundance factor",
					myFormatter.format(proteinTmp.getNSAFNorm()), psiMsCv).getCvParam());
		}
		// distinct peptide sequences
		Set<PSM> psms = DTASelect2MzIdUtil.getPsms(dtaSelectProteins);
		Set<String> sequences = new HashSet<String>();
		for (PSM psm : psms) {
			sequences.add(psm.getSequence());
		}
		int distinctPeptideSequences = sequences.size();
		pdh.getCvParam().add(DTASelect2MzIdUtil.getCVParam(DISTINCT_PEP_SEQUENCSE, "distinct peptide sequences",
				String.valueOf(distinctPeptideSequences), psiMsCv).getCvParam());
		for (Protein protein : dtaSelectProteins) {
			final Set<PSM> psMs = protein.getPSMs();
			for (PSM dtaSelectPSM : psMs) {
				final PeptideHypothesis peptideHypothesis = getPeptideHypothesis(protein, dtaSelectPSM);
				if (!pdh.getPeptideHypothesis().contains(peptideHypothesis)) {
					pdh.getPeptideHypothesis().add(peptideHypothesis);
				}
			}
		}
		return pdh;
	}

	private String getPDHID(DBSequence dbSequence) {
		return dbSequence.getId().replace("DBSeq_", "");
	}

	private PeptideHypothesis getPeptideHypothesis(Protein dtaSelectProtein, PSM dtaSelectPSM) throws IOException {
		final SpectrumIdentificationItemRef siiref = getSpectrumIdentificationItemRef(dtaSelectPSM);
		final PeptideEvidence pe = getPeptideEvidence(dtaSelectPSM, dtaSelectProtein, null);
		if (!phByPeptideEvidence.containsKey(pe.getId())) {
			PeptideHypothesis ph = new PeptideHypothesis();
			ph.setPeptideEvidence(pe);
			phByPeptideEvidence.put(pe.getId(), ph);
		}
		final PeptideHypothesis peptideHypothesis = phByPeptideEvidence.get(pe.getId());
		boolean found = false;
		for (SpectrumIdentificationItemRef siiref2 : peptideHypothesis.getSpectrumIdentificationItemRef()) {
			if (siiref2.getSpectrumIdentificationItemRef().equals(siiref.getSpectrumIdentificationItemRef())) {
				found = true;
			}
		}
		if (!found) {
			peptideHypothesis.getSpectrumIdentificationItemRef().add(siiref);
		}
		return peptideHypothesis;
	}

	private SpectrumIdentificationItemRef getSpectrumIdentificationItemRef(PSM dtaSelectPSM) throws IOException {
		SpectrumIdentificationItemRef siref = new SpectrumIdentificationItemRef();
		siref.setSpectrumIdentificationItem(getSpectrumIdentificationItem(dtaSelectPSM));
		return siref;
	}

	private FragmentationTable getFragmentationTable() {
		// TODO Auto-generated method stub
		return null;
	}

	private Long getNumSeqSearched() {
		// TODO Auto-generated method stub
		return null;
	}

	private Inputs getInputs() throws IOException {
		Inputs ret = new Inputs();
		ret.getSearchDatabase().add(getFastaDBSequence());
		int i = 1;
		for (File file : dtaSelectFiles) {
			ret.getSourceFile().add(getSourceFile("DTASelect_" + i++, file));
		}
		for (LabeledSearchType lst : LabeledSearchType.values()) {
			final List<InputSpectra> inputSpectraList = getInputSpectraList(lst);
			Set<String> spectradataIds = new HashSet<String>();
			for (InputSpectra inputSpectra : inputSpectraList) {
				final SpectraData spectraData = inputSpectra.getSpectraData();
				if (!spectradataIds.contains(spectraData.getId())) {
					ret.getSpectraData().add(spectraData);
					spectradataIds.add(spectraData.getId());
				}
			}
		}
		return ret;
	}

	private SourceFile getSourceFile(String id, File file) {
		SourceFile ret = new SourceFile();
		ret.setFileFormat(getDtaSelectFileFormat());
		ret.setId(id);
		ret.setLocation(file.getAbsolutePath());
		ret.setName(FilenameUtils.getName(file.getAbsolutePath()));
		try {
			String md5Checksum = MD5Checksum.getMD5Checksum(file.getAbsolutePath());
			ret.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1000568", "MD5", md5Checksum, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		} catch (Exception e) {
			e.printStackTrace();
		}

		return ret;
	}

	private AnalysisProtocolCollection getAnalysisProcotolCollection() throws IOException {
		AnalysisProtocolCollection ret = new AnalysisProtocolCollection();
		ret.setProteinDetectionProtocol(getProteinDetectionProtocol());
		for (LabeledSearchType lst : LabeledSearchType.values()) {
			final SpectrumIdentificationProtocol spectrumIdentificationProtocol = getSpectrumIdentificationProtocol(
					lst);
			if (spectrumIdentificationProtocol != null) {
				ret.getSpectrumIdentificationProtocol().add(spectrumIdentificationProtocol);
			}
		}

		return ret;
	}

	private AnalysisCollection getAnalysisCollection() throws IOException {
		AnalysisCollection ret = new AnalysisCollection();
		ret.setProteinDetection(getProteinDetection());
		for (String spectraFileName : dtaSelectParser.getSpectraFileNames()) {
			final SpectrumIdentification spectrumIdentification = getSpectrumIdentification(spectraFileName);
			// check if not already there
			if (!ret.getSpectrumIdentification().contains(spectrumIdentification)) {
				ret.getSpectrumIdentification().add(spectrumIdentification);
			}
		}
		return ret;
	}

	private SpectrumIdentification getSpectrumIdentification(String fileID) throws IOException {
		final LabeledSearchType lst = getLabeledSearchTypeByFileName(fileID);
		final String siID = "SI_" + fileID + "_" + lst.getKey();

		if (!siMap.containsKey(siID)) {
			final SpectrumIdentificationProtocol sip = getSpectrumIdentificationProtocol(lst);
			// only if there is a sip, create the si
			if (sip != null) {
				SpectrumIdentification si = new SpectrumIdentification();
				si.setId(siID);
				final SpectrumIdentificationList sil = getSpectrumIdentificationList(fileID, lst);
				if (sil != null) {
					si.setSpectrumIdentificationList(sil);
				}
				si.setSpectrumIdentificationProtocol(sip);
				si.getInputSpectra().addAll(getInputSpectraList(fileID, lst));
				si.getSearchDatabaseRef().add(getSearchDatabaseRef());
				if (sil != null) {
					siMap.put(siID, si);
				}
			}
		}
		return siMap.get(siID);
	}

	private SearchDatabaseRef getSearchDatabaseRef() throws IOException {
		SearchDatabaseRef ret = new SearchDatabaseRef();
		ret.setSearchDatabase(getFastaDBSequence());
		return ret;
	}

	private List<InputSpectra> getInputSpectraList(String fileID, LabeledSearchType lst) throws IOException {
		List<InputSpectra> ret = new ArrayList<InputSpectra>();
		Set<String> spectraFileNames = dtaSelectParser.getSpectraFileNames();
		Set<String> spectraDataIds = new HashSet<String>();
		for (String spectraFileName : spectraFileNames) {
			if (spectraFileName.equals(fileID)) {
				final LabeledSearchType labeledSearchTypeByFileName = getLabeledSearchTypeByFileName(spectraFileName);
				if (labeledSearchTypeByFileName != LabeledSearchType.LIGHT) {
					// TODO
					// spectraFileName = spectraFileName.substring(1,
					// spectraFileName.length());
				}
				final InputSpectra inputSpetra = getInputSpetra(spectraFileName);
				if (!spectraDataIds.contains(inputSpetra)) {
					ret.add(inputSpetra);
					spectraDataIds.add(inputSpetra.getSpectraDataRef());
				}
			}
		}

		return ret;
	}

	private List<InputSpectra> getInputSpectraList(LabeledSearchType lst) throws IOException {
		List<InputSpectra> ret = new ArrayList<InputSpectra>();
		Set<String> spectraFileNames = dtaSelectParser.getSpectraFileNames();
		Set<String> spectraDataIds = new HashSet<String>();
		for (String spectraFileName2 : spectraFileNames) {
			if (getLabeledSearchTypeByFileName(spectraFileName2) == lst) {
				if (!spectraDataIds.contains(spectraFileName2)) {
					final InputSpectra inputSpetra = getInputSpetra(spectraFileName2);
					ret.add(inputSpetra);
					spectraDataIds.add(inputSpetra.getSpectraDataRef());
				}
			}
		}

		return ret;
	}

	private InputSpectra getInputSpetra(String spectraFileName) throws IOException {
		InputSpectra ret = new InputSpectra();
		ret.setSpectraData(getSpectraData(spectraFileName));
		return ret;
	}

	private SpectrumIdentificationProtocol getSpectrumIdentificationProtocol(LabeledSearchType lst) throws IOException {
		if (!sipByLST.containsKey(lst)) {
			SpectrumIdentificationProtocol sip = new SpectrumIdentificationProtocol();
			sip.setId("SIP_" + lst.getKey());
			sip.setAnalysisSoftware(getSearchEngine());
			sip.setEnzymes(SearchParametersManager.getInstance().getEnzymes(lst));
			sip.setFragmentTolerance(SearchParametersManager.getInstance().getFragmentTolerance(lst));
			sip.setParentTolerance(SearchParametersManager.getInstance().getParentTolerance(lst));
			sip.setSearchType(getMSMSSearchType());
			final ModificationParams modificationParams = SearchParametersManager.getInstance()
					.getModificationParams(lst);
			if (modificationParams != null) {
				sip.setModificationParams(modificationParams);
			}
			final ParamList spectrumThreshold = getSpectrumThreshold();
			if (spectrumThreshold != null) {
				sip.setThreshold(spectrumThreshold);
			} else {
				ParamList paramList = new ParamList();
				paramList.getCvParam().add(DTASelect2MzIdUtil
						.getCVParam("MS:1001494", "no threshold", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
				sip.setThreshold(paramList);
			}
			final ParamList additionalSearchParams = SearchParametersManager.getInstance()
					.getAdditionalSearchParams(lst, getInputSpectraList(lst).size());
			sip.setAdditionalSearchParams(additionalSearchParams);
			// if it is light, it is ok to not having search parameters
			// otherwise, return null
			if (lst != LabeledSearchType.LIGHT && additionalSearchParams == null) {
				return null;
			}

			sip.getMassTable().add(getMassTable(lst));

			// add it at the end, if everything is fine here
			sipByLST.put(lst, sip);
		}
		return sipByLST.get(lst);
	}

	private MassTable getMassTable(LabeledSearchType lst) {
		// create a map with the used modifications
		Map<String, Float> fixedModifications = new HashMap<String, Float>();
		// final ModificationParams modificationParams2 =
		// getModificationParams(lst);
		// for (SearchModification searchModification :
		// modificationParams2.getSearchModification()) {
		// if (searchModification.isFixedMod()) {
		// final List<String> residues = searchModification.getResidues();
		// for (String residue : residues) {
		// fixedModifications.put(residue, searchModification.getMassDelta());
		// }
		// }
		// }

		MassTable ret = new MassTable();
		ret.setId("MT_" + lst.getKey());
		ret.setName("mass table");
		ret.getMsLevel().add(2);
		String aminoacids = "GASPVTCLIXNOBDQKZEMHFRYW";
		AssignMass.getInstance(true);
		for (int index = 0; index < aminoacids.length() - 1; index++) {
			final char aa = aminoacids.charAt(index);
			if (fixedModifications.containsKey(String.valueOf(aa))) {
				ret.getResidue()
						.add(getResidue(aa, AssignMass.getMass(aa) + fixedModifications.get(String.valueOf(aa))));
			} else {
				ret.getResidue().add(getResidue(aa, AssignMass.getMass(aa)));
			}
		}
		return ret;
	}

	private Residue getResidue(char aa, double mass) {
		Residue ret = new Residue();
		ret.setCode(String.valueOf(aa));
		ret.setMass(Double.valueOf(mass).floatValue());
		return ret;
	}

	/**
	 *
	 * @param residues
	 *            null for N-term and C-term modifications
	 * @param location
	 *            0 for N-term mod. 1 for 1st AA, etc...
	 * @param monoMassDelta
	 * @param avgMassDelta
	 * @param fixed
	 * @param nterm
	 * @param cterm
	 * @return
	 */
	private Modification getModification(String residues, Integer location, Double monoMassDelta, Double avgMassDelta,
			boolean nterm, boolean cterm) {
		Modification ret = new Modification();
		ret.setMonoisotopicMassDelta(monoMassDelta);
		ret.setAvgMassDelta(avgMassDelta);
		ret.setLocation(location);
		if (residues != null) {
			for (int i = 0; i < residues.length(); i++) {
				String aa = String.valueOf(residues.charAt(i));
				ret.getResidues().add(aa);
			}
		}
		PeptideModificationUtil peptideModUtil = new PeptideModificationUtil(monoMassDelta, residues);
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

	private ParamList getSpectrumThreshold() throws IOException {

		if (dtaSelectParser.getCommandLineParameter().getParametersMap().containsKey("--fp")) {
			ParamList ret = new ParamList();
			ret.getCvParam()
					.add(DTASelect2MzIdUtil.getCVParam("MS:1001448", "Pep:FDR threshold",
							dtaSelectParser.getCommandLineParameter().getParameterValue("--fp"),
							DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
			return ret;
		}
		return null;
	}

	private File getInputFolder() {
		return new File(dtaSelectFiles.get(0).getParent());
	}

	private Param getMSMSSearchType() {
		ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(MSMSSearchType,
				SearchType.getInstance(cvManager));
		return DTASelect2MzIdUtil.getCVParam(cvTerm, null);
	}

	private ProteinDetection getProteinDetection() throws IOException {
		ProteinDetection ret = new ProteinDetection();
		ret.setId("PD");
		ret.setProteinDetectionList(getProteinDetectionList());
		ret.setProteinDetectionProtocol(getProteinDetectionProtocol());
		for (String fileID : dtaSelectParser.getSpectraFileNames()) {
			for (LabeledSearchType lst : LabeledSearchType.values()) {
				final InputSpectrumIdentifications inputSpectrumIdentification = getInputSpectrumIdentification(fileID,
						lst);
				if (inputSpectrumIdentification != null) {
					ret.getInputSpectrumIdentifications().add(inputSpectrumIdentification);
				}
			}
		}
		return ret;
	}

	private InputSpectrumIdentifications getInputSpectrumIdentification(String fileID, LabeledSearchType lst)
			throws IOException {
		InputSpectrumIdentifications ret = new InputSpectrumIdentifications();
		final SpectrumIdentificationList spectrumIdentificationList = getSpectrumIdentificationList(fileID, lst);
		if (spectrumIdentificationList == null) {
			return null;
		}
		ret.setSpectrumIdentificationList(spectrumIdentificationList);
		return ret;
	}

	private SpectrumIdentificationList getSpectrumIdentificationList(String fileID, LabeledSearchType lst)
			throws IOException {
		final String silID = "SIL_" + fileID + "_" + getSearchEngine().getId() + "_" + lst.getKey();
		if (!spectrumIdentificationListMap.containsKey(silID)) {
			SpectrumIdentificationList spectrumIdentificationList = new SpectrumIdentificationList();
			spectrumIdentificationList.setId(silID);

			spectrumIdentificationList.setNumSequencesSearched(getNumSeqSearched());

			final Map<String, PSM> dtaSelectPSMsByPSMID = dtaSelectParser.getPSMsByPSMID();
			boolean atLeastOnePSM = false;
			for (PSM psm : dtaSelectPSMsByPSMID.values()) {
				// include the PSM only if it was searched in the
				// prefractionation step
				if (psm.getMSRun().getRunId().equals(fileID)) {
					// include the PSM only if it was searched in the
					// corresponding
					// labeledSearchType
					if (getLabeledSearchTypeByFileName(psm.getMSRun().getRunId()) == lst) {

						spectrumIdentificationList.getSpectrumIdentificationResult()
								.add(getSpectrumIdentificationResult(psm));
						atLeastOnePSM = true;
					}
				} else {
					log.debug("ignoring psm in " + lst.name() + "_ " + lst.getKey());
				}
			}

			if (!atLeastOnePSM) {
				return null;
			}

			// tag as final PSM list
			// removed since at https://github.com/HUPO-PSI/mzIdentML/issues/5
			// HUPO-PSI decided that only final results are going to be encoded
			// spectrumIdentificationList.getCvParam()
			// .add(DTASelect2MzIdUtil
			// .getCVParam(FINAL_PSM_LIST_CV, "final PSM list", null,
			// DTASelect2MzIdUtil.getPSIMsCv())
			// .getCvParam());

			spectrumIdentificationList.setFragmentationTable(getFragmentationTable());

			log.info("SpectrumIdentificationList associated with pre-fractionation step '" + fileID + "' and search '"
					+ lst.name() + "' created.");

			// add to the map at the end
			spectrumIdentificationListMap.put(silID, spectrumIdentificationList);

		}
		return spectrumIdentificationListMap.get(silID);
	}

	private SpectrumIdentificationResult getSpectrumIdentificationResult(PSM dtaSelectPSM) throws IOException {

		final String sir_id = getSpectrumIdentificationResultID(dtaSelectPSM);
		if (sirMap.containsKey(sir_id)) {
			return sirMap.get(sir_id);
		} else {
			SpectrumIdentificationResult sir = new SpectrumIdentificationResult();
			sirMap.put(sir_id, sir);

			sir.setId(sir_id);
			sir.setSpectrumID(getSpectrumID(dtaSelectPSM));

			final String runId = dtaSelectPSM.getMSRun().getRunId();
			final SpectraData spectraData = getSpectraData(runId);
			sir.setSpectraData(spectraData);
			sir.getSpectrumIdentificationItem().add(getSpectrumIdentificationItem(dtaSelectPSM));
			// spectrum title
			sir.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1000796", "spectrum title",
					dtaSelectPSM.getPSMIdentifier(), DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
			return sir;
		}
	}

	private String getSpectrumIdentificationResultID(PSM dtaSelectPSM) {
		final String spectrumID = getSpectrumID(dtaSelectPSM);
		String spectraFileName = dtaSelectPSM.getMSRun().getRunId();
		final String sir_id = spectraFileName + "_Spec_" + spectrumID;
		return sir_id;
	}

	private String getSpectrumID(PSM dtaSelectPSM) {
		String spectraFileName = dtaSelectPSM.getMSRun().getRunId();
		if (getLabeledSearchTypeByFileName(spectraFileName) != LabeledSearchType.LIGHT) {
			spectraFileName = spectraFileName.substring(1, spectraFileName.length());
		}
		String spectrumID = null;
		if (ms2ReaderByFileName.containsKey(spectraFileName)) {
			final MS2Reader ms2Reader = ms2ReaderByFileName.get(spectraFileName);

			if (referenceToSpectra == ReferenceToSpectra.MS2) {
				final String key = dtaSelectPSM.getScanNumber() + "." + dtaSelectPSM.getScanNumber() + "."
						+ dtaSelectPSM.getChargeState();
				final Integer spectrumIndex = ms2Reader.getSpectrumIndexByScan(key);
				if (spectrumIndex != null) {
					spectrumID = "index=" + (spectrumIndex);
				}
			} else if (referenceToSpectra == ReferenceToSpectra.MZXML) {
				spectrumID = "scan=" + dtaSelectPSM.getScanNumber();
			}
			if (spectrumID != null) {
				return spectrumID;
			}

			log.warn("Spectrum not found in ms2  '" + ms2Reader.getFileName() + "' for PSM '"
					+ dtaSelectPSM.getPSMIdentifier() + "' scan='" + dtaSelectPSM.getScanNumber() + "' charge='"
					+ dtaSelectPSM.getChargeState() + "'");

		} else {
			if (referenceToSpectra == ReferenceToSpectra.MZXML) {
				spectrumID = "scan=" + dtaSelectPSM.getScanNumber();
			} else {
				if (!spectraFileNamesNotFound.contains(spectraFileName) && !ignoreSpectra
						&& referenceToSpectra == ReferenceToSpectra.MS2) {
					log.warn("MS2 file reader not found for '" + spectraFileName + "' file reference");
					spectraFileNamesNotFound.add(spectraFileName);
				}
			}
			if (spectrumID != null) {
				return spectrumID;
			}
		}
		return dtaSelectPSM.getPSMIdentifier();
	}

	public LabeledSearchType getLabeledSearchTypeByFileName(String fileName) {
		if (labeledSearchTypeByFileName.containsKey(fileName)) {
			return labeledSearchTypeByFileName.get(fileName);
		}
		String path = getDTASelectFolder(dtaSelectFiles);

		// substract the first letter
		String lightfile = path + File.separator + fileName.substring(1, fileName.length()) + ".ms2";
		File lf = new File(lightfile);
		if (lf.exists()) {
			// depending on the first letter will be one labeled search type or
			// the other
			String key = fileName.substring(0, 1);
			for (LabeledSearchType labeledSearchType : LabeledSearchType.values()) {
				if (labeledSearchType.getKey().equals(key)) {
					labeledSearchTypeByFileName.put(fileName, labeledSearchType);
					return labeledSearchType;
				}
			}
		}
		// by default
		labeledSearchTypeByFileName.put(fileName, LabeledSearchType.LIGHT);
		return LabeledSearchType.LIGHT;

	}

	private Double getScanStartTime(PSM dtaSelectPSM) {
		final String spectraFileName = dtaSelectPSM.getMSRun().getRunId();
		if (ms2ReaderByFileName.containsKey(spectraFileName)) {
			final MS2Reader ms2Reader = ms2ReaderByFileName.get(spectraFileName);
			final String key = dtaSelectPSM.getScanNumber() + "." + dtaSelectPSM.getScanNumber() + "."
					+ dtaSelectPSM.getChargeState();
			final Double rt = ms2Reader.getSpectrumRTByScan(key);
			if (rt != null) {
				return rt;
			}
		}
		return null;
	}

	private SpectrumIdentificationItem getSpectrumIdentificationItem(PSM dtaSelectPSM) throws IOException {
		final String siiID = getSpectrumIdentificationResultID(dtaSelectPSM) + "_"
				+ getPeptide(dtaSelectPSM, null).getId();
		if (siiMap.containsKey(siiID)) {
			final SpectrumIdentificationItem spectrumIdentificationItem = siiMap.get(siiID);
			return spectrumIdentificationItem;
		} else {
			SpectrumIdentificationItem ret = new SpectrumIdentificationItem();
			ret.setId(siiID);
			siiMap.put(siiID, ret);

			final Integer charge = Integer.valueOf(dtaSelectPSM.getChargeState());
			final double nonMonoExperimentalMZ = MassesUtil.getMassToCharge(dtaSelectPSM.getExperimentalMH(), charge);

			final double calculatedMZ = MassesUtil.getMassToCharge(dtaSelectPSM.getCalcMH(), charge);
			ret.setCalculatedMassToCharge(calculatedMZ);

			// final Double monoMassExperimental =
			// getMonoMassFromNonMonoMass(nonMonoExperimentalMZ, calculatedMZ,
			// charge);
			ret.setExperimentalMassToCharge(nonMonoExperimentalMZ);

			try {
				ret.setChargeState(charge);
			} catch (NumberFormatException e) {

			}
			if (dtaSelectPSM.getPI() != null) {
				ret.setCalculatedPI(dtaSelectPSM.getPI().floatValue());
			}

			ret.setPassThreshold(true);
			ret.setPeptide(getPeptide(dtaSelectPSM, null));
			ret.setRank(1);
			ret.setSample(getSample(false));

			// peptide evidences
			final Set<Protein> proteins = dtaSelectPSM.getProteins();
			Set<String> peptideEvidenceIds = new HashSet<String>();
			for (Protein dtaSelectProtein : proteins) {
				final PeptideEvidenceRef peptideEvidenceRef = getPeptideEvidenceRef(dtaSelectPSM, dtaSelectProtein);
				if (!peptideEvidenceIds.contains(peptideEvidenceRef.getPeptideEvidenceRef())) {
					ret.getPeptideEvidenceRef().add(peptideEvidenceRef);
					peptideEvidenceIds.add(peptideEvidenceRef.getPeptideEvidenceRef());
				}
			}
			// scores
			ret.getCvParam().add(getDeltaCnScore(getDeltaCn(dtaSelectPSM)));
			ret.getCvParam().add(getXCorrScore(getXCorr(dtaSelectPSM)));

			// TODO
			// THIS HAS BEEN ADDED IN ORDER TO BE COMPATIBLE WITH SKYLINE UNTIL
			// THEY MODIFY THEIR SOFTWARE IN ORDER TO INCLUDE EVERYTHING IN THE
			// FILE
			if (skylineCompatible) {
				ret.getCvParam().add(getFAKE_SCORE_FOR_SKYLINE());
			}

			// TODO check if valid
			final CvParam scanStartTimeCV = getScanStartTimeCV(dtaSelectPSM);
			if (scanStartTimeCV != null) {
				ret.getCvParam().add(scanStartTimeCV);
			}
			return ret;
		}
	}

	private Double getMonoMassFromNonMonoMass(Double nonMonoMZ, Double monoMZToApproach, Integer charge) {
		if (nonMonoMZ == null || monoMZToApproach == null || charge == null) {
			return null;
		}
		double diff = Math.abs(nonMonoMZ - monoMZToApproach);
		double monoMZApproached = nonMonoMZ;
		while (true) {
			monoMZApproached -= AssignMass.DIFFMASSC12C13 / charge;
			double tmpDiff = Math.abs(monoMZApproached - monoMZToApproach);
			if (tmpDiff > diff) {
				final double ret = monoMZApproached + AssignMass.DIFFMASSC12C13 / charge;
				return ret;
			}
			diff = tmpDiff;
		}
	}

	private Double getNonMonoMassFromMonoMass(Double nonMono, Double mono, Integer charge) {
		if (nonMono == null || mono == null || charge == null) {
			return null;
		}
		double diff = Math.abs(nonMono - mono);
		double nonMonoTmp = mono;
		while (true) {
			nonMonoTmp += AssignMass.DIFFMASSC12C13 / charge;
			double tmpDiff = Math.abs(nonMonoTmp - nonMono);
			if (tmpDiff > diff) {
				final double ret = nonMonoTmp - AssignMass.DIFFMASSC12C13 / charge;
				return ret;
			}
			diff = tmpDiff;
		}
	}

	private CvParam getScanStartTimeCV(PSM psm) {
		Double scanStartTime = getScanStartTime(psm);
		if (scanStartTime != null) {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(SCAN_START_TIME_CV,
					PeptideResultsDetailsCVSet.getInstance(cvManager));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, myFormatter.format(scanStartTime));
			return cvParam.getCvParam();
		}
		return null;
	}

	private CvParam getFAKE_SCORE_FOR_SKYLINE() {
		final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(FAKE_SCORE_FOR_SKYLINE,
				Score.getInstance(cvManager));
		// give allways the value 1
		final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, "100");
		return cvParam.getCvParam();
	}

	private Double getDeltaCn(PSM psm) {
		try {
			if (psm != null) {
				final Set<edu.scripps.yates.utilities.proteomicsmodel.Score> scores = psm.getScores();
				for (edu.scripps.yates.utilities.proteomicsmodel.Score score : scores) {
					if (score.getScoreName().toLowerCase().contains("deltacn")) {
						return Double.valueOf(score.getValue());
					}
				}
			}
		} catch (NumberFormatException e) {
		}
		return null;
	}

	private Double getXCorr(PSM psm) {
		try {
			if (psm != null) {
				final Set<edu.scripps.yates.utilities.proteomicsmodel.Score> scores = psm.getScores();
				for (edu.scripps.yates.utilities.proteomicsmodel.Score score : scores) {
					if (score.getScoreName().toLowerCase().contains("xcorr")) {
						return Double.valueOf(score.getValue());
					}
				}
			}
		} catch (NumberFormatException e) {
		}
		return null;
	}

	private CvParam getDeltaCnScore(Double deltacn) throws IOException {
		if (deltacn == null) {
			return null;
		}

		if (dtaSelectParser.getSearchEngines().contains(DTASelectParser.SEQUEST)) {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(DELTACN_SEQUEST,
					Score.getInstance(cvManager));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, String.valueOf(deltacn));
			return cvParam.getCvParam();
		} else if (dtaSelectParser.getSearchEngines().contains(DTASelectParser.PROLUCID)) {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(DELTACN_PROLUCID,
					Score.getInstance(cvManager));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, String.valueOf(deltacn));
			return cvParam.getCvParam();
		}
		return null;
	}

	private CvParam getXCorrScore(Double xcorr) throws IOException {
		if (xcorr == null) {
			return null;
		}
		final ControlVocabularyManager om = DTASelect2MzIdUtil.getOntologyManager();

		if (dtaSelectParser.getSearchEngines().contains(DTASelectParser.SEQUEST)) {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(XCORR_SEQUEST,
					Score.getInstance(om));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, String.valueOf(xcorr));
			return cvParam.getCvParam();
		} else if (dtaSelectParser.getSearchEngines().contains(DTASelectParser.PROLUCID)) {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(XCORR_PROLUCID,
					Score.getInstance(om));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, String.valueOf(xcorr));
			return cvParam.getCvParam();
		}
		return null;
	}

	private PeptideEvidenceRef getPeptideEvidenceRef(PSM dtaSelectPSM, Protein dtaSelectProtein) throws IOException {
		PeptideEvidenceRef ret = new PeptideEvidenceRef();
		ret.setPeptideEvidence(getPeptideEvidence(dtaSelectPSM, dtaSelectProtein, null));
		return ret;
	}

	private Sample getSample(boolean createSampleIfNoExists) {
		if (sample == null && createSampleIfNoExists) {
			sample = new Sample();
			sample.setId("Sample");
			sample.setName("Sample name");
			sample.getContactRole().add(getYatesLabContactRole(DTASelect2MzIdUtil.getResearcherRole()));
		}
		return sample;
	}

	private SpectraData getSpectraData(String fileID) throws IOException {
		final LabeledSearchType labeledSearchTypeByFileName = getLabeledSearchTypeByFileName(fileID);
		if (labeledSearchTypeByFileName != LabeledSearchType.LIGHT) {
			// TODO
			// fileID = fileID.substring(1, fileID.length());
		}
		if (spectraDataBySpectraFileName.containsKey(fileID)) {
			return spectraDataBySpectraFileName.get(fileID);
		}
		SpectraData spectraData = new SpectraData();
		spectraDataBySpectraFileName.put(fileID, spectraData);
		String extension = null;
		if (referenceToSpectra == ReferenceToSpectra.MS2) {
			spectraData.setFileFormat(getMS2FileFormat());
			extension = ".ms2";
			spectraData.setExternalFormatDocumentation("http://fields.scripps.edu/sequest/SQTFormat.html");
		} else if (referenceToSpectra == ReferenceToSpectra.MZXML) {
			spectraData.setFileFormat(getMzXMLFileFormat());
			extension = ".mzXML";
		}
		spectraData.setId(fileID);
		final String filePath = dtaSelectParser.getRunPath() + File.separator + fileID + extension;
		spectraData.setLocation(filePath);
		spectraData.setSpectrumIDFormat(getSpectrumIDFormat());

		return spectraData;
	}

	private SpectrumIDFormat getSpectrumIDFormat() {
		SpectrumIDFormat ret = new SpectrumIDFormat();
		if (referenceToSpectra == ReferenceToSpectra.MS2) {
			ret.setCvParam(DTASelect2MzIdUtil.getCVParam("MS:1000774", "multiple peak list nativeID format", null,
					DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		} else if (referenceToSpectra == ReferenceToSpectra.MZXML) {
			ret.setCvParam(DTASelect2MzIdUtil
					.getCVParam("MS:1000776", "scan number only nativeID format", null, DTASelect2MzIdUtil.getPSIMsCv())
					.getCvParam());
		}
		// Used for conversion of peak list files with multiple spectra, i.e.
		// MGF, PKL, merged DTA files. Index is the spectrum number in the file,
		// starting from 0.
		// Native format defined by index=xsd:nonNegativeInteger.
		return ret;
	}

	private FileFormat getMS2FileFormat() {
		FileFormat format = new FileFormat();
		format.setCvParam(DTASelect2MzIdUtil
				.getCVParam(MS2_FORMAT_CV, "MS2 format", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		return format;
	}

	private FileFormat getMzXMLFileFormat() {
		FileFormat format = new FileFormat();
		format.setCvParam(DTASelect2MzIdUtil
				.getCVParam(MzXML_FORMAT_CV, "ISB mzXML format", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		return format;
	}

	private FileFormat getDtaSelectFileFormat() {
		FileFormat format = new FileFormat();
		format.setCvParam(DTASelect2MzIdUtil
				.getCVParam(DTASELECT_FORMAT_CV, "DTASelect format", null, DTASelect2MzIdUtil.getPSIMsCv())
				.getCvParam());
		return format;
	}

	private ProteinDetectionProtocol getProteinDetectionProtocol() throws IOException {
		if (pdp == null) {
			pdp = new ProteinDetectionProtocol();
			pdp.setAnalysisSoftware(getSearchEngine());
			pdp.setId("PDP");
			final ParamList proteinThreshold = getProteinThreshold();
			if (proteinThreshold != null) {
				pdp.setThreshold(proteinThreshold);
			} else {
				ParamList paramList = new ParamList();
				paramList.getCvParam().add(DTASelect2MzIdUtil
						.getCVParam("MS:1001494", "no threshold", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
				pdp.setThreshold(paramList);
			}
		}
		return pdp;
	}

	private ParamList getProteinThreshold() throws IOException {
		if (dtaSelectParser.getCommandLineParameter().getParametersMap().containsKey("--fp")) {
			ParamList ret = new ParamList();
			ret.getCvParam()
					.add(DTASelect2MzIdUtil.getCVParam("MS:1001447", "prot:FDR threshold",
							dtaSelectParser.getCommandLineParameter().getParameterValue("--fp"),
							DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
			return ret;
		}
		return null;
	}

	private ProteinDetectionList getProteinDetectionList() throws IOException {
		if (pdl == null) {
			pdl = new ProteinDetectionList();
			pdl.setId("PDL_" + getSearchEngine().getId());
			final List<ProteinGroup> proteinGroups = getGroups();
			int numPAG = 0;
			for (ProteinGroup proteinGroup : proteinGroups) {
				pdl.getProteinAmbiguityGroup().add(getProteinAmbiguityGroup(proteinGroup));
				numPAG++;
			}
			// number of proteins

			pdl.getCvParam()
					.add(DTASelect2MzIdUtil.getCVParam(COUNT_IDENTIFIED_PROTEINS_CV, "count of identified proteins",
							String.valueOf(numPAG), DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		}
		return pdl;
	}

	private SequenceCollection getSequenceCollection() throws IOException {
		SequenceCollection ret = new SequenceCollection();
		final HashMap<String, Set<Protein>> proteinMap = dtaSelectParser.getProteins();
		if (proteinMap != null) {
			for (Set<Protein> proteinSet : proteinMap.values()) {
				for (Protein protein : proteinSet) {
					final DBSequence dbSequence = getDBSequence(protein);
					// check if already there
					if (!ret.getDBSequence().contains(dbSequence)) {
						ret.getDBSequence().add(dbSequence);
					}
					final Set<PSM> psMs = protein.getPSMs();
					for (PSM dtaSelectPSM : psMs) {
						PeptideEvidence peptideEvidence = getPeptideEvidence(dtaSelectPSM, protein, ret);
						if (!ret.getPeptideEvidence().contains(peptideEvidence)) {
							ret.getPeptideEvidence().add(peptideEvidence);
						}
						final Peptide peptide = getPeptide(dtaSelectPSM, ret);
						if (!ret.getPeptide().contains(peptide)) {
							ret.getPeptide().add(peptide);
						}
					}

				}
			}

		}
		return ret;
	}

	private PeptideEvidence getPeptideEvidence(PSM dtaSelectPSM, Protein protein, SequenceCollection sequenceCollection)
			throws IOException {
		String key = getPeptideEvidenceKey(dtaSelectPSM, protein);
		if (peptideEvidencesByKey.containsKey(key)) {
			return peptideEvidencesByKey.get(key);
		}
		PeptideEvidence peptideEvidence = new PeptideEvidence();
		peptideEvidencesByKey.put(key, peptideEvidence);
		if (sequenceCollection != null) {
			sequenceCollection.getPeptideEvidence().add(peptideEvidence);
		}

		final DBSequence dbSequence = getDBSequence(protein);
		peptideEvidence.setDBSequence(dbSequence);
		Peptide peptide = getPeptide(dtaSelectPSM, sequenceCollection);
		peptideEvidence.setId(key);
		peptideEvidence.setPeptide(peptide);
		peptideEvidence.setPre(String.valueOf(dtaSelectPSM.getBeforeSeq()));
		peptideEvidence.setPost(String.valueOf(dtaSelectPSM.getAfterSeq()));
		// check if decoy
		final Pattern pattern = getDecoyPattern();
		if (pattern != null) {
			final Matcher matcher = pattern.matcher(dbSequence.getAccession());
			if (matcher.find()) {
				peptideEvidence.setIsDecoy(true);
			} else {
				peptideEvidence.setIsDecoy(false);
			}
		}
		return peptideEvidence;
	}

	private Pattern getDecoyPattern() throws IOException {
		if (decoyPattern == null) {
			final String decoyPattern2 = getDecoyPatternString();
			if (decoyPattern2 != null) {
				decoyPattern = Pattern.compile(decoyPattern2);
			}
		}
		return decoyPattern;
	}

	private String getPeptideEvidenceKey(PSM dtaSelectPSM, Protein protein) throws IOException {
		return "Pep_" + getPeptide(dtaSelectPSM, null).getId() + "_DBSeq" + getDBSequence(protein).getId();
	}

	private Peptide getPeptide(PSM dtaSelectPSM, SequenceCollection sequenceCollection) throws IOException {

		// the peptide KEY is the sequence + the modifications

		Peptide peptide = new Peptide();

		peptide.setPeptideSequence(dtaSelectPSM.getSequence());
		// diff modifications
		if (dtaSelectPSM.getPTMs() != null && !dtaSelectPSM.getPTMs().isEmpty()) {
			for (PTM modification : dtaSelectPSM.getPTMs()) {
				for (PTMSite site : modification.getPTMSites()) {
					peptide.getModification().add(getModification(modification, site));
				}
			}
		}
		// fix modifications. They are treated as different masses in the AA
		// in the search of PROLucid
		List<Modification> fixedModifications = getFixedModifications(dtaSelectPSM);
		peptide.getModification().addAll(fixedModifications);
		String peptideID = getPeptideKey(peptide);
		peptide.setId(peptideID);

		if (!peptides.containsKey(peptideID)) {
			peptides.put(peptideID, peptide);
		} else {
			peptide = peptides.get(peptideID);
		}
		if (sequenceCollection != null && !sequenceCollection.getPeptide().contains(peptide)) {
			sequenceCollection.getPeptide().add(peptide);
		}
		return peptide;
	}

	private String getPeptideKey(Peptide peptide) {
		StringBuilder sb = new StringBuilder();
		final String peptideSequence = peptide.getPeptideSequence();
		final List<Modification> modifications = peptide.getModification();
		Collections.sort(modifications, new Comparator<Modification>() {

			@Override
			public int compare(Modification o1, Modification o2) {
				int location1 = o1.getLocation() != null ? o1.getLocation() : -1;
				int location2 = o2.getLocation() != null ? o2.getLocation() : -1;
				return Integer.compare(location1, location2);
			}
		});
		// n term is location=0
		for (Modification modification : modifications) {
			if (modification.getLocation() != null && modification.getLocation() == 0) {
				sb.append("(" + myFormatter4digits.format(modification.getMonoisotopicMassDelta()) + ")");
			}
		}
		for (int index = 0; index < peptideSequence.length(); index++) {
			final int position = index + 1;
			sb.append(peptideSequence.charAt(index));
			for (Modification modification : modifications) {
				if (modification.getLocation() != null && modification.getLocation() == position) {
					sb.append("(" + myFormatter4digits.format(modification.getMonoisotopicMassDelta()) + ")");
				}
			}
		}
		// c term is location=peptide length+1
		for (Modification modification : modifications) {
			if (modification.getLocation() != null && modification.getLocation() == peptideSequence.length() + 1) {
				sb.append("(" + myFormatter4digits.format(modification.getMonoisotopicMassDelta()) + ")");
			}
		}
		return sb.toString();
	}

	private List<PTM> getSortedByName(Collection<PTM> ptms) {
		List<PTM> list = new ArrayList<PTM>();
		list.addAll(ptms);
		Collections.sort(list, new Comparator<PTM>() {

			@Override
			public int compare(PTM o1, PTM o2) {
				return o1.getName().compareTo(o2.getName());
			}
		});
		return list;
	}

	private List<Modification> getFixedModifications(PSM dtaSelectPSM) throws IOException {
		List<Modification> ret = new ArrayList<Modification>();
		if (dtaSelectPSM != null) {
			LabeledSearchType lst = getLabeledSearchTypeByFileName(dtaSelectPSM.getMSRun().getRunId());
			final SearchXmlFile searchParameters = SearchParametersManager.getInstance().getSearchParameters(lst);
			if (searchParameters != null) {

				// N term fix
				if (searchParameters.getNTermStaticMod() != null) {
					final double massDiff = Double.valueOf(searchParameters.getNTermStaticMod());
					if (Double.compare(massDiff, 0.0) != 0) {
						ret.add(getModification(null, null, massDiff, null, true, false));
					}
				}

				// C term fix
				if (searchParameters.getCTermStaticMod() != null) {
					final double massDiff = Double.valueOf(searchParameters.getCTermStaticMod());
					if (Double.compare(massDiff, 0.0) != 0) {
						ret.add(getModification(null, null, massDiff, null, false, true));
					}
				}

				// modifications fix
				if (searchParameters.getStaticmods() != null) {

					// loop over the peptide sequence
					final String sequence = dtaSelectPSM.getSequence();
					for (int index = 0; index < sequence.length(); index++) {
						final char aa = sequence.charAt(index);
						for (String diffAndResidue : searchParameters.getStaticmods()) {

							double massDiff = Double.valueOf(diffAndResidue.split(" ")[0]);
							String residues = diffAndResidue.split(" ")[1];
							for (int index2 = 0; index2 < residues.length(); index2++) {
								final char aa2 = residues.charAt(index2);
								if (aa == aa2) {
									ret.add(getModification(residues, index + 1, massDiff, null, false, false));
								}
							}

						}
					}
				}

			}
		}
		return ret;
	}

	private Modification getModification(PTM modification, PTMSite site) {
		if (modification == null || site == null) {
			return null;
		}
		PeptideModificationUtil peptideModUtil = new PeptideModificationUtil(modification, site);
		Modification ret = new Modification();
		ret.setLocation(peptideModUtil.getPosition());
		ret.setMonoisotopicMassDelta(peptideModUtil.getMonoDelta());
		ret.getResidues().add(peptideModUtil.getResidues());
		if (peptideModUtil.getAccession() != null) {
			ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam(peptideModUtil.getAccession(), peptideModUtil.getName(),
					null, DTASelect2MzIdUtil.getUnimodCv()).getCvParam());
		} else {
			ret.getCvParam()
					.add(DTASelect2MzIdUtil
							.getCVParam("MS:1001460", "Unknown modification", null, DTASelect2MzIdUtil.getPSIMsCv())
							.getCvParam());
		}

		return ret;

	}

	private DBSequence getDBSequence(Protein dtaSelectProtein) throws IOException {
		// look into the map
		if (dbSequences.containsKey(dtaSelectProtein.getAccession())) {
			return dbSequences.get(dtaSelectProtein.getAccession());
		}
		DBSequence dbSequence = new DBSequence();
		dbSequences.put(dtaSelectProtein.getAccession(), dbSequence);
		final Pair<String, String> acc = FastaParser.getACC(dtaSelectProtein.getAccession());
		dbSequence.setAccession(acc.getFirstelement());
		dbSequence.setId(getDBSequenceID(acc.getFirstelement()));
		dbSequence.setSearchDatabase(getFastaDBSequence());
		final String description = FastaParser.getDescription(dtaSelectProtein.getPrimaryAccession().getDescription());
		dbSequence.setName(description);
		dbSequence.getCvParam().add(DTASelect2MzIdUtil
				.getCVParam(PROTEIN_DESCRIPTION_CV, "protein description", description, DTASelect2MzIdUtil.getPSIMsCv())
				.getCvParam());
		return dbSequence;
	}

	private String getDBSequenceID(String acc) {
		return "DBSeq_" + acc;
	}

	private SearchDatabase getFastaDBSequence() throws IOException {
		if (searchDatabase == null) {
			searchDatabase = new SearchDatabase();
			String fastaPath = dtaSelectParser.getFastaPath();
			searchDatabase.setId("FASTA");
			searchDatabase.setName(FilenameUtils.getBaseName(fastaPath));
			searchDatabase.setLocation(fastaPath);
			searchDatabase.setDatabaseName(DTASelect2MzIdUtil.getUserParam(FilenameUtils.getBaseName(fastaPath), null));
			final Cv psiMsCv = DTASelect2MzIdUtil.getPSIMsCv();

			searchDatabase.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1001015", "database original uri", new File(fastaPath).toURI().toString(), psiMsCv)
					.getCvParam());
			if (isDecoy()) {
				searchDatabase.getCvParam().add(DTASelect2MzIdUtil
						.getCVParam("MS:1001197", "DB composition target+decoy", null, psiMsCv).getCvParam());
				searchDatabase.getCvParam()
						.add(DTASelect2MzIdUtil
								.getCVParam("MS:1001283", "decoy DB accession regexp", getDecoyPatternString(), psiMsCv)
								.getCvParam());
				if (getDecoyPatternString().equalsIgnoreCase("reverse")) {
					searchDatabase.getCvParam().add(DTASelect2MzIdUtil
							.getCVParam("MS:1001195", "decoy DB type reverse", null, psiMsCv).getCvParam());
				} else if (getDecoyPatternString().toLowerCase().contains("rnd")
						|| getDecoyPatternString().toLowerCase().contains("random")) {
					searchDatabase.getCvParam().add(DTASelect2MzIdUtil
							.getCVParam("MS:1001196", "decoy DB type randomized", null, psiMsCv).getCvParam());
				} else if (getDecoyPatternString().toLowerCase().contains("shuffle")) {
					searchDatabase.getCvParam().add(DTASelect2MzIdUtil
							.getCVParam("MS:1001452", "decoy DB type shuffle", null, psiMsCv).getCvParam());
				}
			}
			searchDatabase.setFileFormat(getFastaFormat());
		}
		return searchDatabase;
	}

	private boolean isDecoy() throws IOException {
		if (getDecoyPatternString() != null)
			return true;
		return false;
	}

	private String getDecoyPatternString() throws IOException {
		final String decoyPattern2 = dtaSelectParser.getDecoyPattern();
		if (decoyPattern2 == null) {
			return ".*Reverse.*";
		}
		return decoyPattern2;
	}

	private FileFormat getFastaFormat() {
		FileFormat ret = new FileFormat();
		ret.setCvParam(DTASelect2MzIdUtil
				.getCVParam("MS:1001348", "FASTA format", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		return ret;
	}

	private AnalysisSampleCollection getAnalysisSampleCollection(boolean createSampleIfNoExists) {
		AnalysisSampleCollection ret = new AnalysisSampleCollection();
		final Sample sample2 = getSample(createSampleIfNoExists);
		if (sample2 == null) {
			return null;
		}
		ret.getSample().add(sample2);
		return ret;
	}

	private AuditCollection getAuditCollection() {
		AuditCollection ret = new AuditCollection();
		ret.getOrganization().add(getYatesLabContact());
		ret.getPersonOrOrganization().add(getScrippsOrganization().getOrganization());
		return ret;
	}

	private Provider getDefaultProvider() {
		Provider ret = new Provider();
		ret.setContactRole(getYatesLabContactRole(DTASelect2MzIdUtil.getResearcherRole()));
		ret.setId("PROVID");
		return ret;
	}

	private AnalysisSoftwareList getAnalysisSoftwareList() throws IOException {
		AnalysisSoftwareList ret = new AnalysisSoftwareList();

		ret.getAnalysisSoftware().add(getSearchEngine());
		ret.getAnalysisSoftware().add(getDTASelectSoftware());
		return ret;
	}

	private AnalysisSoftware getDTASelectSoftware() throws IOException {
		AnalysisSoftware ret = new AnalysisSoftware();
		ret.setId("DTASelect");
		ret.setName("DTASelect analysis software");
		ret.setVersion(dtaSelectParser.getDTASelectVersion());
		ret.setUri("http://fields.scripps.edu/DTASelect/");
		ret.setContactRole(
				getYatesLabContactRole(DTASelect2MzIdUtil.getSoftwareVendorRole("Yates lab software developers")));
		// TODO change by DTASELECT when supported by skyline
		if (skylineCompatible) {
			final ControlVocabularyTerm dtaSelectCV = DTASelect2MzIdUtil.getCvTermByAcc(SEQUEST_CV,
					SoftwareName.getInstance(DTASelect2MzIdUtil.getOntologyManager()));
			ret.setSoftwareName(DTASelect2MzIdUtil.getCVParam(dtaSelectCV, null));
		} else {
			ret.setSoftwareName(DTASelect2MzIdUtil.getCVParam(DTASELECT_TOOL_CV, "DTASelect", null,
					DTASelect2MzIdUtil.getPSIMsCv()));
		}

		return ret;
	}

	private AnalysisSoftware getSearchEngine() throws IOException {
		AnalysisSoftware ret = null;
		if (isProLuCIDSearch()) {
			return getProLuCIDAnalysisSoftware();
		}
		if (isSequestSearch()) {
			return getSequestAnalysisSoftware();
		}
		return ret;
	}

	private AnalysisSoftware getSequestAnalysisSoftware() throws IOException {
		AnalysisSoftware ret = new AnalysisSoftware();
		ret.setId(DTASelectParser.SEQUEST);
		ret.setName("SEQUEST search engine");
		ret.setVersion(dtaSelectParser.getSearchEngineVersion());
		ret.setUri("http://fields.scripps.edu/sequest/");
		ret.setContactRole(
				getYatesLabContactRole(DTASelect2MzIdUtil.getSoftwareVendorRole("Yates lab software developers")));
		final ControlVocabularyTerm sequestCV = DTASelect2MzIdUtil.getCvTermByAcc(SEQUEST_CV,
				SoftwareName.getInstance(DTASelect2MzIdUtil.getOntologyManager()));
		ret.setSoftwareName(DTASelect2MzIdUtil.getCVParam(sequestCV, null));
		return ret;
	}

	private AnalysisSoftware getProLuCIDAnalysisSoftware() throws IOException {
		AnalysisSoftware ret = new AnalysisSoftware();
		ret.setId(DTASelectParser.PROLUCID);
		ret.setName("ProLuCID search engine");
		ret.setVersion(dtaSelectParser.getSearchEngineVersion());
		ret.setUri("http://fields.scripps.edu/downloadfile2.php?name=ProLuCID&filename=&id=12");
		ret.setContactRole(
				getYatesLabContactRole(DTASelect2MzIdUtil.getSoftwareVendorRole("Yates lab software developers")));
		// TODO change by sequest
		if (skylineCompatible) {
			final ControlVocabularyTerm prolucidCV = DTASelect2MzIdUtil.getCvTermByAcc(SEQUEST_CV,
					SoftwareName.getInstance(DTASelect2MzIdUtil.getOntologyManager()));
			ret.setSoftwareName(DTASelect2MzIdUtil.getCVParam(prolucidCV, null));
		} else {
			final Param prolucidCV = DTASelect2MzIdUtil.getCVParam(PROLUCID_CV, "ProLuCID", null,
					DTASelect2MzIdUtil.getPSIMsCv());
			ret.setSoftwareName(prolucidCV);
		}
		return ret;
	}

	private ContactRole getYatesLabContactRole(Role role) {
		ContactRole ret = new ContactRole();
		ret.setContact(getYatesLabContact());
		ret.setRole(role);
		return ret;
	}

	private Organization getYatesLabContact() {
		Organization ret = new Organization();
		ret.setId(YATES_LAB_CONTACT_ID);
		ret.setName("Yates lab");
		// TODO add organization cv terms
		ret.setParent(getScrippsOrganization());
		ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1000587", "contact address",
				"10550 North Torrey Pines Road, SR302", DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		ret.getCvParam().add(DTASelect2MzIdUtil
				.getCVParam("MS:1000589", "contact email", "salvador@scripps.edu", DTASelect2MzIdUtil.getPSIMsCv())
				.getCvParam());
		return ret;
	}

	private ParentOrganization getScrippsOrganization() {
		ParentOrganization ret = new ParentOrganization();
		Organization scrippsOrg = new Organization();
		scrippsOrg.setId(SCRIPPS_ORG_ID);
		scrippsOrg.setName("The Scripps Research Institute");
		// scrippsOrg.getCvParam().add(DTASelect2MzIdUtil
		// .getCVParam("MS:1001755", "contact phone number", "858-784-1000",
		// DTASelect2MzIdUtil.getPSIMsCv())
		// .getCvParam());
		scrippsOrg.getCvParam()
				.add(DTASelect2MzIdUtil
						.getCVParam("MS:1000587", "contact address",
								"10550 North Torrey Pines Road, La Jolla, CA 92037", DTASelect2MzIdUtil.getPSIMsCv())
						.getCvParam());
		ret.setOrganization(scrippsOrg);
		return ret;
	}

	private boolean isProLuCIDSearch() throws IOException {
		return isSearchEngine(DTASelectParser.PROLUCID);
	}

	private boolean isSequestSearch() throws IOException {
		return isSearchEngine(DTASelectParser.SEQUEST);
	}

	private boolean isSearchEngine(String searchEngineName) throws IOException {
		if (dtaSelectParser != null) {
			final Set<String> searchEngines = dtaSelectParser.getSearchEngines();
			for (String searchEngine : searchEngines) {
				if (searchEngine.equals(searchEngineName)) {
					return true;
				}
			}
		}
		return false;
	}

	public static void main(String[] args) {

		setupCommandLineOptions();
		File inputFile = null;
		String decoyRegexp = null;
		boolean uniqueOutput = false;
		String inputFileName = "DTASelect-filter.txt";
		boolean recursiveInputFileSearch = false;
		boolean ignoreSpectra = false;
		MzIdentMLVersion version = MzIdentMLVersion.VERSION_1_1;
		boolean skylineCompatible = true;
		ReferenceToSpectra referenceToSpectra = ReferenceToSpectra.MS2;
		CommandLineParser parser = new BasicParser();
		try {
			CommandLine cmd = parser.parse(options, args);
			if (cmd.hasOption("i")) {
				inputFile = new File(cmd.getOptionValue("i"));
			} else {
				errorInParameters("Input file is missing");
			}

			if (cmd.hasOption("d")) {
				decoyRegexp = cmd.getOptionValue("d");
				try {
					Pattern.compile(decoyRegexp);
				} catch (PatternSyntaxException e) {
					errorInParameters("Decoy regular expression '" + decoyRegexp + "' is not well formed.");
				}
				log.info("DECOY regular expression provided: '" + decoyRegexp
						+ "'. DECOY hits will be ignored in the output.");
			} else {
				log.info(
						"No DECOY regular expression provided. If DECOY matches are included in the DTASelect file, they would also be included in the output mzIdentML file");
			}
			if (cmd.hasOption("ns")) {
				ignoreSpectra = true;
				log.info("Ignoring MS2 files. No spectra will be readed.");
			}

			if (cmd.hasOption("rs")) {
				try {
					referenceToSpectra = ReferenceToSpectra.valueOf(cmd.getOptionValue("rs"));

				} catch (Exception e) {
					errorInParameters(
							"Non valid value for 'rs' option. Valid values are: " + ReferenceToSpectra.getCSVString());
				}
			}
			if (!ignoreSpectra) {
				log.info("Reference to spectra in " + referenceToSpectra.name() + " format");
			}
			if (cmd.hasOption("r")) {
				recursiveInputFileSearch = true;
				if (inputFile.isDirectory()) {
					log.info("Option 'r' detected. Look for input files recursively");
				} else {
					errorInParameters("Option 'r' can only be used with a folder as input in option 'i'");
				}
			} else {
				log.info("No Option 'r' detected. Only one input file will be converted.");
			}
			if (cmd.hasOption("n")) {
				inputFileName = cmd.getOptionValue("n");
				if (inputFile.isFile()) {
					errorInParameters("Option 'n' can only be used with a folder as input in option 'i'");
				}
			}

			if (cmd.hasOption("u")) {
				uniqueOutput = true;
				log.info("Option 'u' = '" + uniqueOutput + "' detected. A unique mzIdentML file will be created");
			}
			if (uniqueOutput) {
				log.info("Creating a single mzIdentML file from all input files.");
			} else {
				log.info("Creating one mzIdentML file from each input file.");
			}

			if (cmd.hasOption("v")) {
				String v = cmd.getOptionValue("v");
				final MzIdentMLVersion versionByString = MzIdentMLVersion.getByString(v);
				if (versionByString == null) {
					errorInParameters("Version '" + v + "' is not supported. Versions supported are: "
							+ MzIdentMLVersion.getCSVString());
				} else {
					version = versionByString;
				}
				log.info("Option 'v' = '" + v + "' detected. Converting to mzIdentML version '"
						+ version.getVersionString() + "'");
			} else {
				log.info("Option 'v' not detected. Converting by default to mzIdentML version '"
						+ version.getVersionString() + "'");
			}
			if (cmd.hasOption("sky")) {
				final String sky = cmd.getOptionValue("sky");
				try {
					skylineCompatible = Boolean.valueOf(sky);
				} catch (Exception e) {
					errorInParameters("'" + sky + "' is not a valid boolean value");
				}
				log.info("Option 'sky' = '" + sky + "' detected. Compatibility with skyline is '" + skylineCompatible
						+ "'");
			} else {
				log.info("Option 'sky' not detected. Compatibility with skyline is by default '" + skylineCompatible
						+ "'");
			}

			List<File> inputFiles = getInputFiles(inputFile, inputFileName, recursiveInputFileSearch);
			if (uniqueOutput) {
				DTASelect2MzId conversor = new DTASelect2MzId(inputFiles, new File(System.getProperty("user.dir")
						+ File.separator + FilenameUtils.getBaseName(inputFileName) + ".mzid"), version);
				conversor.setDecoyRegexp(decoyRegexp);
				conversor.setIgnoreSpectra(ignoreSpectra);
				conversor.setSkylineCompatible(skylineCompatible);
				conversor.setReferenceToSpectra(referenceToSpectra);
				conversor.convert();
			} else {
				for (File file : inputFiles) {
					inputFileName = FilenameUtils.getBaseName(file.getAbsolutePath());
					final File outputMzIdentMLFile = new File(
							file.getParentFile().getAbsolutePath() + File.separator + inputFileName + ".mzid");
					log.info("Using file name: '" + inputFileName + "' for output file(s).");
					DTASelect2MzId conversor = new DTASelect2MzId(file, outputMzIdentMLFile, version);
					conversor.setDecoyRegexp(decoyRegexp);
					conversor.setIgnoreSpectra(ignoreSpectra);
					conversor.setSkylineCompatible(skylineCompatible);
					conversor.setReferenceToSpectra(referenceToSpectra);

					conversor.convert();
				}
			}
		} catch (ParseException e) {
			errorInParameters(e.getMessage());
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println(e.getMessage());
		}

	}

	private void setSkylineCompatible(boolean skylineCompatible) {
		this.skylineCompatible = skylineCompatible;
	}

	/**
	 * @return the skylineCompatible
	 */
	public boolean isSkylineCompatible() {
		return skylineCompatible;
	}

	private void setIgnoreSpectra(boolean ignoreSpectra) {
		this.ignoreSpectra = ignoreSpectra;
	}

	private static List<File> getInputFiles(File inputFile, String inputFileName, boolean recursiveInputFileSearch) {
		List<File> ret = new ArrayList<File>();
		if (inputFile.isFile()) {
			ret.add(inputFile);
		} else {
			if (!recursiveInputFileSearch) {
				File file = new File(inputFile.getAbsolutePath() + File.separator + inputFileName);
				ret.add(file);
			} else {
				List<File> files = findFiles(inputFile, inputFileName);
				ret.addAll(files);
			}
		}
		String plural = ret.size() > 1 ? "s" : "";
		log.info(ret.size() + " file" + plural + " found:");
		int i = 1;
		for (File file : ret) {
			log.info(i++ + "- " + file.getAbsolutePath());
		}
		return ret;
	}

	private static List<File> findFiles(File inputFolder, String inputFileName) {
		List<File> ret = new ArrayList<File>();
		// look in current folder
		File file = new File(inputFolder.getAbsolutePath() + File.separator + inputFileName);
		if (file.exists()) {
			ret.add(file);
		}
		// look in subfolders
		final File[] listFiles = inputFolder.listFiles();
		for (File file2 : listFiles) {
			if (file2.isDirectory()) {
				ret.addAll(findFiles(file2, inputFileName));
			}
		}
		return ret;
	}

	private void setDecoyRegexp(String decoyRegexp) {
		decoyPatternString = decoyRegexp;
	}

	private static void errorInParameters(String header) {
		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		if (header == null) {
			formatter.printHelp("DTASelect2MzId -i [input file or folder]", options);
		} else {
			formatter.printHelp("DTASelect2MzId -i [input file or folder]",
					"\n\n************\n" + header + "\n************\n\n", options,
					"Contact Salvador Martinez-Bartolome at salvador@scripps.edu for more help");
		}
		System.exit(0);
	}

	private static void setupCommandLineOptions() {
		// create Options object
		options = new Options();

		// add t option
		options.addOption("i", "input", true, "path to the input file (or folder)");
		options.addOption("d", "decoy", true,
				"[OPTIONAL] decoy regular expression. Ignores matching entries. Example: 'Reverse'.");
		options.addOption("r", "recursive", false,
				"[OPTIONAL] In case of using a folder as input '-i', it will search recursively for all the DTASelect-filter.txt files.");
		options.addOption("n", "file_name", true,
				"[OPTIONAL] To use a input file name different from the default 'DTASelect-filter.txt'");
		options.addOption("u", "unique_output_file", false,
				"[OPTIONAL] A single mzIdentML file will be created collapsing all the information of all the input files."
						+ " Otherwise, a mzIdentML file will be created for each input file.");
		options.addOption("ns", "no_spectra", false,
				"[OPTIONAL] If provided, no MS2 files will be readed in order to match PSMs with spectra in the output file");
		options.addOption("rs", "referenceToSpectra", true, "[OPTIONAL] Reference to spectra. Possible values: "
				+ ReferenceToSpectra.getCSVString() + ". Default: '" + ReferenceToSpectra.MS2 + "'");
		options.addOption("v", "version", true, "[OPTIONAL] Version of the output mzIdentML '-v', Possible values: "
				+ MzIdentMLVersion.getCSVString() + ". Default: " + MzIdentMLVersion.VERSION_1_1.getVersionString());
		options.addOption("sky", "skyline", true,
				"[OPTIONAL] Whether the generated mzIdentML will be compatible with skyline for importing results. If 'true', "
						+ "the file will contain additional features to make it compatible. Possible values: 'true' or 'false'. Default: 'true'");
	}

}
