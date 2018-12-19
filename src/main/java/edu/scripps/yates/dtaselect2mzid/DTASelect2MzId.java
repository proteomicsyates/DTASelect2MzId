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

import edu.scripps.yates.dtaselect2mzid.util.DTASelect2MzIdUtil;
import edu.scripps.yates.dtaselect2mzid.util.LabeledSearchType;
import edu.scripps.yates.dtaselect2mzid.util.MS2Reader;
import edu.scripps.yates.dtaselect2mzid.util.MzIdentMLVersion;
import edu.scripps.yates.dtaselect2mzid.util.PeptideModificationUtil;
import edu.scripps.yates.dtaselect2mzid.util.ReferenceToSpectra;
import edu.scripps.yates.dtaselect2mzid.util.SpectrumIdentificationResultDetailsCVSet;
import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.utilities.checksum.MD5Checksum;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.masses.AssignMass;
import edu.scripps.yates.utilities.masses.MassesUtil;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.factories.ProteinEx;
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
	private DTASelectParser dtaSelectParser;
	private final AssignMass assignMass = AssignMass.getInstance(true);
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
		for (final File file : dtaSelectFiles) {
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

		dtaSelectParser = new DTASelectParser(dtaSelectFiles);
		dtaSelectParser.setDecoyPattern(decoyPatternString);

		if (!ignoreSpectra && referenceToSpectra == ReferenceToSpectra.MS2) {
			final String path = getDTASelectFolder(dtaSelectFiles);
			final Set<String> spectraFileFullPaths = dtaSelectParser.getSpectraFileFullPaths();
			for (final String spectraFileFullPath : spectraFileFullPaths) {
				// MS2Reader
				try {
					File ms2File = new File(spectraFileFullPath);
					if (!ms2File.exists()) {
						ms2File = new File(path + File.separator + spectraFileFullPath + ".ms2");
					}
					if (ms2File.exists()) {
						final MS2Reader ms2Reader = new MS2Reader(ms2File);
						log.info("Setting up ms2 reader for file: '" + ms2File.getAbsolutePath());
						String spectraFileName = FilenameUtils.getBaseName(ms2File.getAbsolutePath());
						if (ms2ReaderByFileName.containsKey(spectraFileName)) {
							log.info(ms2ReaderByFileName.get(spectraFileName).getFileName());
						}

						if (getLabeledSearchTypeByFileName(spectraFileName) != LabeledSearchType.LIGHT) {
							spectraFileName = spectraFileName.substring(1, spectraFileName.length());
						}
						ms2ReaderByFileName.put(spectraFileName, ms2Reader);
					}
				} catch (final IllegalArgumentException e) {
					log.warn(e.getMessage());
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
		MzIdentMLMarshaller m = null;
		try {
			writer = new FileWriter(output);
			m = new MzIdentMLMarshaller();
			// XML header
			log.info("Creating XML header...");
			writer.write(m.createXmlHeader() + "\n");
			// mzIdentML start tag
			writer.write(createMzIdentMLStartTag(m, mzIdentMLVersion, "DTASelect2MzId") + "\n");

			log.info("Creating CvList element...");
			final CvList cvList = DTASelect2MzIdUtil.getCVList();
			m.marshal(cvList, writer);
			writer.write("\n");

			//
			log.info("Creating AnalysisSoftwareList element...");
			final AnalysisSoftwareList analysisSoftwareList = getAnalysisSoftwareList();
			m.marshal(analysisSoftwareList, writer);
			writer.write("\n");

			log.info("Creating Provider element...");
			final Provider provider = getDefaultProvider();
			m.marshal(provider, writer);
			writer.write("\n");

			//
			//
			log.info("Creating AuditCollection element...");
			final AuditCollection auditCollection = getAuditCollection();
			m.marshal(auditCollection, writer);
			writer.write("\n");
			//
			log.info("Creating AnalysisSampleCollection element...");
			final AnalysisSampleCollection analysisSampleCollection = getAnalysisSampleCollection(true);
			if (analysisSampleCollection != null) {
				m.marshal(analysisSampleCollection, writer);
				writer.write("\n");
			}
			//
			log.info("Creating SequenceCollection element...");
			final SequenceCollection sequenceCollection = getSequenceCollection();
			log.info("Sequence collection with " + sequenceCollection.getDBSequence().size() + " DBSequences, "
					+ sequenceCollection.getPeptide().size() + " Peptides and "
					+ sequenceCollection.getPeptideEvidence().size() + " PeptideEvidences");
			m.marshal(sequenceCollection, writer);
			writer.write("\n");

			log.info("Creating SpectrumIdentificationResults...");
			// get fisrt the spectra
			for (final PSM psm : dtaSelectParser.getPSMsByPSMID().values()) {
				// System.out.println(psm.getMSRun().getRunId() + "\t" +
				// psm.getPSMIdentifier());
				final SpectrumIdentificationResult specIdentRes = getSpectrumIdentificationResult(psm);
			}
			log.info(dtaSelectParser.getPSMsByPSMID().size() + " SpectrumIdentificationResult elements created");

			//
			log.info("Creating AnalysisCollection element...");
			final AnalysisCollection analysisCollection = getAnalysisCollection();
			m.marshal(analysisCollection, writer);
			writer.write("\n");
			//
			log.info("Creating AnalysisProcoloCollection element...");
			final AnalysisProtocolCollection analysisProtocolCollection = getAnalysisProcotolCollection();
			m.marshal(analysisProtocolCollection, writer);
			writer.write("\n");
			writer.flush();
			//
			log.info("Creating DataCollection element...");
			writer.write(m.createDataCollectionStartTag() + "\n");
			//
			log.info("Creating Inputs element...");
			final Inputs inputs = getInputs();
			m.marshal(inputs, writer);
			writer.write("\n");
			//
			log.info("Creating AnalysisData element...");
			writer.write(m.createAnalysisDataStartTag() + "\n");
			//
			// one spectrumIdentificationList by each fileName
			// (pre-fractionation) and each labelledSearchTag if
			// available
			for (final String fileID : dtaSelectParser.getSpectraFileNames()) {
				for (final LabeledSearchType lst : LabeledSearchType.values()) {
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

		} finally {
			writer.write(m.createMzIdentMLClosingTag());
			log.info("File created at: " + output.getAbsolutePath());

			final Set<String> errorMessages = PeptideModificationUtil.getErrorMessages();
			if (!errorMessages.isEmpty()) {
				System.err.println(
						"Some modifications were not identified properly and were translated as 'unknown' modifications:");
				for (final String errorMessage : errorMessages) {
					System.err.println(errorMessage);
				}
			}
			if (writer != null)
				writer.close();
			if (mzIdentMLVersion == MzIdentMLVersion.VERSION_1_2) {
				fixLines();
			}
		}
	}

	private void fixLines() {
		try {
			final Path path = Paths.get(output.toURI());
			final Charset charset = StandardCharsets.UTF_8;

			String content = new String(Files.readAllBytes(path), charset);
			content = content.replaceAll("http://psidev.info/psi/pi/mzIdentML/1.1",
					"http://psidev.info/psi/pi/mzIdentML/1.2");
			Files.write(path, content.getBytes(charset));

		} catch (final IOException e) {
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
		final StringBuffer sb = new StringBuffer();
		final String MZIDML_NS = "http://psidev.info/psi/pi/mzIdentML/1.2";
		final String MZIDML_VERSION = "1.2.0";
		final String MZIDML_SCHEMA = "https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0-candidate.xsd";
		final String MZIDSCHEMA_LOCATION = "http://psidev.info/psi/pi/mzIdentML/1.2 " + MZIDML_SCHEMA;
		// tag opening plus id attribute
		sb.append("<MzIdentML id=\"").append(id).append("\"");
		// further attributes
		sb.append(" version=\"").append(MZIDML_VERSION).append("\"");
		sb.append(" xmlns=\"").append(MZIDML_NS).append("\"");
		sb.append(" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"");
		sb.append(" xsi:schemaLocation=\"").append(MZIDSCHEMA_LOCATION).append("\"");
		final DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss");
		sb.append(" creationDate=\"").append(dfm.format(getCreationDate())).append("\"");
		// finally close the tag
		sb.append(" >");

		return sb.toString();
	}

	private Date getCreationDate() {
		if (dtaSelectFiles != null && !dtaSelectFiles.isEmpty()) {
			for (final File file : dtaSelectFiles) {
				try {
					final BasicFileAttributes attr = Files.readAttributes(file.toPath(), BasicFileAttributes.class);
					final Date date = new Date(attr.creationTime().toMillis());
					return date;
				} catch (final IOException e) {
				}
			}
		}
		return new Date();

	}

	private String getDTASelectFolder(List<File> dtaSelectFiles2) {

		return dtaSelectFiles2.get(0).getParentFile().getAbsolutePath();
	}

	private List<ProteinGroup> getGroups() throws IOException {
		if (groups == null) {
			log.info("Grouping  " + dtaSelectParser.getProteins().size()
					+ " proteins according to PAnalyzer algorithm...");
			groups = new ArrayList<ProteinGroup>();
			final Set<GroupableProtein> groupableProteins = new HashSet<GroupableProtein>();
			for (final Protein protein : dtaSelectParser.getProteins()) {

				groupableProteins.add(protein);
			}
			final PAnalyzer panalyzer = new PAnalyzer(true);
			groups = panalyzer.run(groupableProteins);
			log.info("Proteins grouped in " + groups.size() + " groups.");
		}
		return groups;

	}

	private ProteinAmbiguityGroup getProteinAmbiguityGroup(ProteinGroup group) throws IOException {
		final String pagKey = getProteinAmbiguityGroupKey(group);
		if (pag.containsKey(pagKey)) {
			return pag.get(pagKey);
		}
		final ProteinAmbiguityGroup ret = new ProteinAmbiguityGroup();
		pag.put(pagKey, ret);
		ret.setId(pagKey);
		final Map<String, Set<Protein>> proteinsByAcc = getProteinsByAcc(group);
		final Map<String, ProteinEvidence> evidences = new HashMap<String, ProteinEvidence>();
		for (final String acc : proteinsByAcc.keySet()) {
			evidences.put(acc, proteinsByAcc.get(acc).iterator().next().getEvidence());
		}
		for (final String acc : proteinsByAcc.keySet()) {
			final Set<Protein> proteinSet = proteinsByAcc.get(acc);
			ret.getProteinDetectionHypothesis().add(getProteinDetectionHypothesis(proteinSet));
		}

		ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam(PG_PASS_THRESHOLD_CV, "protein group passes threshold",
				"true", DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		return ret;
	}

	private Map<String, Set<Protein>> getProteinsByAcc(ProteinGroup group) {
		final Map<String, Set<Protein>> map = new HashMap<String, Set<Protein>>();
		for (final GroupableProtein groupableProtein : group) {
			final String accession = groupableProtein.getAccession();
			if (map.containsKey(accession)) {
				map.get(accession).add((Protein) groupableProtein);
			} else {
				final Set<Protein> set = new HashSet<Protein>();
				set.add((Protein) groupableProtein);
				map.put(accession, set);
			}
		}
		return map;
	}

	private String getProteinAmbiguityGroupKey(ProteinGroup group) throws IOException {
		final StringBuilder sb = new StringBuilder();
		final List<String> pdhIds = new ArrayList<String>();
		final Map<String, Set<Protein>> proteinsByAcc = getProteinsByAcc(group);
		for (final Set<Protein> proteinSet : proteinsByAcc.values()) {
			final String pdhId = getProteinDetectionHypothesis(proteinSet).getId();
			if (!pdhIds.contains(pdhId)) {
				pdhIds.add(pdhId);
			}
		}
		for (final String pdhId : pdhIds) {
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
		final Protein protein1 = dtaSelectProteins.iterator().next();
		final DBSequence dbSequence = getDBSequence(protein1);
		final String pdhid = getPDHID(dbSequence);
		if (pdhMapByID.containsKey(pdhid)) {
			return pdhMapByID.get(pdhid);
		}
		final ProteinDetectionHypothesis pdh = new ProteinDetectionHypothesis();
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

		// protein coverage
		pdh.getCvParam().add(DTASelect2MzIdUtil.getCVParam(PROTEIN_COVERAGE_CV, "sequence coverage",
				myFormatter.format(protein1.getCoverage()), psiMsCv).getCvParam());
		// empai value
		if (protein1.getEmpai() != null) {
			pdh.getCvParam()
					.add(DTASelect2MzIdUtil
							.getCVParam(EMPAI_VALUE_CV, "emPAI value", myFormatter.format(protein1.getEmpai()), psiMsCv)
							.getCvParam());
		}
		// NSAF
		pdh.getCvParam().add(DTASelect2MzIdUtil.getCVParam(NSAF_CV, "normalized spectral abundance factor",
				myFormatter.format(protein1.getNsaf_norm()), psiMsCv).getCvParam());

		// distinct peptide sequences
		final Set<PSM> psms = DTASelect2MzIdUtil.getPsms(dtaSelectProteins);
		final Set<String> sequences = new HashSet<String>();
		for (final PSM psm : psms) {
			sequences.add(psm.getSequence());
		}
		final int distinctPeptideSequences = sequences.size();
		pdh.getCvParam().add(DTASelect2MzIdUtil.getCVParam(DISTINCT_PEP_SEQUENCSE, "distinct peptide sequences",
				String.valueOf(distinctPeptideSequences), psiMsCv).getCvParam());
		for (final Protein protein : dtaSelectProteins) {
			final List<PSM> psMs = protein.getPSMs();
			for (final PSM dtaSelectPSM : psMs) {
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
			final PeptideHypothesis ph = new PeptideHypothesis();
			ph.setPeptideEvidence(pe);
			phByPeptideEvidence.put(pe.getId(), ph);
		}
		final PeptideHypothesis peptideHypothesis = phByPeptideEvidence.get(pe.getId());
		boolean found = false;
		for (final SpectrumIdentificationItemRef siiref2 : peptideHypothesis.getSpectrumIdentificationItemRef()) {
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
		final SpectrumIdentificationItemRef siref = new SpectrumIdentificationItemRef();
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
		final Inputs ret = new Inputs();
		ret.getSearchDatabase().add(getFastaDBSequence());
		int i = 1;
		for (final File file : dtaSelectFiles) {
			ret.getSourceFile().add(getSourceFile("DTASelect_" + i++, file));
		}
		for (final LabeledSearchType lst : LabeledSearchType.values()) {
			final List<InputSpectra> inputSpectraList = getInputSpectraList(lst);
			final Set<String> spectradataIds = new HashSet<String>();
			for (final InputSpectra inputSpectra : inputSpectraList) {
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
		final SourceFile ret = new SourceFile();
		ret.setFileFormat(getDtaSelectFileFormat());
		ret.setId(id);
		ret.setLocation(file.getAbsolutePath());
		ret.setName(FilenameUtils.getName(file.getAbsolutePath()));
		try {
			final String md5Checksum = MD5Checksum.getMD5ChecksumFromFileName(file);
			ret.getCvParam().add(DTASelect2MzIdUtil
					.getCVParam("MS:1000568", "MD5", md5Checksum, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		} catch (final Exception e) {
			e.printStackTrace();
		}

		return ret;
	}

	private AnalysisProtocolCollection getAnalysisProcotolCollection() throws IOException {
		final AnalysisProtocolCollection ret = new AnalysisProtocolCollection();
		ret.setProteinDetectionProtocol(getProteinDetectionProtocol());
		for (final LabeledSearchType lst : LabeledSearchType.values()) {
			final SpectrumIdentificationProtocol spectrumIdentificationProtocol = getSpectrumIdentificationProtocol(
					lst);
			if (spectrumIdentificationProtocol != null) {
				ret.getSpectrumIdentificationProtocol().add(spectrumIdentificationProtocol);
			}
		}

		return ret;
	}

	private AnalysisCollection getAnalysisCollection() throws IOException {
		final AnalysisCollection ret = new AnalysisCollection();
		ret.setProteinDetection(getProteinDetection());
		for (final String spectraFileName : dtaSelectParser.getSpectraFileNames()) {
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
				final SpectrumIdentification si = new SpectrumIdentification();
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
		final SearchDatabaseRef ret = new SearchDatabaseRef();
		ret.setSearchDatabase(getFastaDBSequence());
		return ret;
	}

	private List<InputSpectra> getInputSpectraList(String fileID, LabeledSearchType lst) throws IOException {
		final List<InputSpectra> ret = new ArrayList<InputSpectra>();
		final Set<String> spectraFileNames = dtaSelectParser.getSpectraFileNames();
		final Set<String> spectraDataIds = new HashSet<String>();
		for (final String spectraFileName : spectraFileNames) {
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
		final List<InputSpectra> ret = new ArrayList<InputSpectra>();
		final Set<String> spectraFileNames = dtaSelectParser.getSpectraFileNames();
		final Set<String> spectraDataIds = new HashSet<String>();
		for (final String spectraFileName2 : spectraFileNames) {
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
		final InputSpectra ret = new InputSpectra();
		ret.setSpectraData(getSpectraData(spectraFileName));
		return ret;
	}

	private SpectrumIdentificationProtocol getSpectrumIdentificationProtocol(LabeledSearchType lst) throws IOException {
		if (!sipByLST.containsKey(lst)) {
			final SpectrumIdentificationProtocol sip = new SpectrumIdentificationProtocol();
			// set analysis software to DTASelect
			// otherwise PRIDE Inspector fails
			String key = lst.getKey();
			if (key != null && !"".equals(key)) {
				key = "_" + key;
			}
			sip.setId("SIP" + key);
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
				final ParamList paramList = new ParamList();
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
		final Map<String, Float> fixedModifications = new HashMap<String, Float>();
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

		final MassTable ret = new MassTable();
		ret.setId("MT_" + lst.getKey());
		ret.setName("mass table");
		ret.getMsLevel().add(2);
		final String aminoacids = "GASPVTCLIXNOBDQKZEMHFRYW";
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
		final Residue ret = new Residue();
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
		final Modification ret = new Modification();
		ret.setMonoisotopicMassDelta(monoMassDelta);
		ret.setAvgMassDelta(avgMassDelta);
		ret.setLocation(location);
		if (residues != null) {
			for (int i = 0; i < residues.length(); i++) {
				final String aa = String.valueOf(residues.charAt(i));
				ret.getResidues().add(aa);
			}
		}
		final PeptideModificationUtil peptideModUtil = new PeptideModificationUtil(monoMassDelta, residues);
		if (peptideModUtil.getAccession() != null) {
			ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam(peptideModUtil.getAccession(), peptideModUtil.getName(),
					null, DTASelect2MzIdUtil.getPsimodCv()).getCvParam());
		} else {
			ret.getCvParam()
					.add(DTASelect2MzIdUtil
							.getCVParam("MS:1001460", "unknown modification", null, DTASelect2MzIdUtil.getPSIMsCv())
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
			final ParamList ret = new ParamList();
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
		final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(MSMSSearchType,
				SearchType.getInstance(cvManager));
		return DTASelect2MzIdUtil.getCVParam(cvTerm, null);
	}

	private ProteinDetection getProteinDetection() throws IOException {
		final ProteinDetection ret = new ProteinDetection();
		ret.setId("PD");
		ret.setProteinDetectionList(getProteinDetectionList());
		ret.setProteinDetectionProtocol(getProteinDetectionProtocol());
		for (final String fileID : dtaSelectParser.getSpectraFileNames()) {
			for (final LabeledSearchType lst : LabeledSearchType.values()) {
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
		final InputSpectrumIdentifications ret = new InputSpectrumIdentifications();
		final SpectrumIdentificationList spectrumIdentificationList = getSpectrumIdentificationList(fileID, lst);
		if (spectrumIdentificationList == null) {
			return null;
		}
		ret.setSpectrumIdentificationList(spectrumIdentificationList);
		return ret;
	}

	private String getNotNullUnderscoreAndSearchEngine() throws IOException {
		final AnalysisSoftware searchEngine = getSearchEngine();
		if (searchEngine != null && searchEngine.getId() != null) {
			return "_" + searchEngine.getId();
		}
		return "";
	}

	private SpectrumIdentificationList getSpectrumIdentificationList(String fileID, LabeledSearchType lst)
			throws IOException {
		String silID = "SIL_" + fileID + getNotNullUnderscoreAndSearchEngine();
		if (!"".equals(lst.getKey())) {
			silID += "_" + lst.getKey();
		}
		if (!spectrumIdentificationListMap.containsKey(silID)) {
			final SpectrumIdentificationList spectrumIdentificationList = new SpectrumIdentificationList();
			spectrumIdentificationList.setId(silID);

			spectrumIdentificationList.setNumSequencesSearched(getNumSeqSearched());

			final Map<String, PSM> dtaSelectPSMsByPSMID = dtaSelectParser.getPSMsByPSMID();
			boolean atLeastOnePSM = false;
			for (final PSM psm : dtaSelectPSMsByPSMID.values()) {
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

			log.info("SpectrumIdentificationList '" + silID + "' created.");

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
			final SpectrumIdentificationResult sir = new SpectrumIdentificationResult();
			sirMap.put(sir_id, sir);

			sir.setId(sir_id);
			sir.setSpectrumID(getSpectrumID(dtaSelectPSM));

			final String runId = dtaSelectPSM.getMSRun().getRunId();
			final SpectraData spectraData = getSpectraData(runId);
			sir.setSpectraData(spectraData);
			sir.getSpectrumIdentificationItem().add(getSpectrumIdentificationItem(dtaSelectPSM));
			// spectrum title
			sir.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1000796", "spectrum title",
					dtaSelectPSM.getIdentifier(), DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
			return sir;
		}
	}

	private String getSpectrumIdentificationResultID(PSM dtaSelectPSM) {
		final String spectrumID = getSpectrumID(dtaSelectPSM);
		final String spectraFileName = dtaSelectPSM.getMSRun().getRunId();
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
				final String key = Integer.valueOf(dtaSelectPSM.getScanNumber()) + "."
						+ Integer.valueOf(dtaSelectPSM.getScanNumber()) + "."
						+ Integer.valueOf(dtaSelectPSM.getChargeState());
				final Integer spectrumIndex = ms2Reader.getSpectrumIndexByScan(key);
				if (spectrumIndex != null) {
					spectrumID = "index=" + (spectrumIndex);
				} else {
					log.warn("Spectrum not found in ms2  '" + ms2Reader.getFileName() + "'. Scan number='"
							+ dtaSelectPSM.getScanNumber() + "' charge state='" + dtaSelectPSM.getChargeState() + "'");
				}
			} else if (referenceToSpectra == ReferenceToSpectra.MZXML) {
				spectrumID = "scan=" + dtaSelectPSM.getScanNumber();
			}
			if (spectrumID != null) {
				return spectrumID;
			}

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
		return dtaSelectPSM.getIdentifier();
	}

	public LabeledSearchType getLabeledSearchTypeByFileName(String fileName) {
		if (labeledSearchTypeByFileName.containsKey(fileName)) {
			return labeledSearchTypeByFileName.get(fileName);
		}
		final String path = getDTASelectFolder(dtaSelectFiles);

		// substract the first letter
		final String lightfile = path + File.separator + fileName.substring(1, fileName.length()) + "."
		// + referenceToSpectra.name().toLowerCase();\
		// modification on Jun 21 2018
		// after talking with Robin we just decided to only check for the
		// presence of ms2 files with the "H", because if someone uses mzXML,
		// then there is no option for Heavy searches
				+ "ms2";
		final File lf = new File(lightfile);
		if (lf.exists()) {
			// depending on the first letter will be one labeled search type or
			// the other
			final String key = fileName.substring(0, 1);
			for (final LabeledSearchType labeledSearchType : LabeledSearchType.values()) {
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
			final SpectrumIdentificationItem ret = new SpectrumIdentificationItem();
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
			} catch (final NumberFormatException e) {

			}
			if (dtaSelectPSM.getPi() != null) {
				ret.setCalculatedPI(dtaSelectPSM.getPi().floatValue());
			}

			ret.setPassThreshold(true);
			ret.setPeptide(getPeptide(dtaSelectPSM, null));
			ret.setRank(1);
			ret.setSample(getSample(false));

			// peptide evidences
			final Set<Protein> proteins = dtaSelectPSM.getProteins();
			final Set<String> peptideEvidenceIds = new HashSet<String>();
			for (final Protein dtaSelectProtein : proteins) {
				final PeptideEvidenceRef peptideEvidenceRef = getPeptideEvidenceRef(dtaSelectPSM, dtaSelectProtein);
				if (!peptideEvidenceIds.contains(peptideEvidenceRef.getPeptideEvidenceRef())) {
					ret.getPeptideEvidenceRef().add(peptideEvidenceRef);
					peptideEvidenceIds.add(peptideEvidenceRef.getPeptideEvidenceRef());
				}
			}
			// scores
			final Double deltaCn = getDeltaCn(dtaSelectPSM);
			if (deltaCn != null) {
				ret.getCvParam().add(getDeltaCnScore(deltaCn));
			}

			final Double xCorr = getXCorr(dtaSelectPSM);
			if (xCorr != null) {
				ret.getCvParam().add(getXCorrScore(xCorr));
			}
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
			final double tmpDiff = Math.abs(monoMZApproached - monoMZToApproach);
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
			final double tmpDiff = Math.abs(nonMonoTmp - nonMono);
			if (tmpDiff > diff) {
				final double ret = nonMonoTmp - AssignMass.DIFFMASSC12C13 / charge;
				return ret;
			}
			diff = tmpDiff;
		}
	}

	private CvParam getScanStartTimeCV(PSM psm) {
		final Double scanStartTime = getScanStartTime(psm);
		if (scanStartTime != null) {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(SCAN_START_TIME_CV,
					SpectrumIdentificationResultDetailsCVSet.getInstance(cvManager));
			if (cvTerm != null) {
				final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, myFormatter.format(scanStartTime));
				return cvParam.getCvParam();
			} else {
				log.warn("scan_start_time CV is not found");
			}
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
				for (final edu.scripps.yates.utilities.proteomicsmodel.Score score : scores) {
					if (score.getScoreName().toLowerCase().contains("deltacn")) {
						return Double.valueOf(score.getValue());
					}
				}
			}
		} catch (final NumberFormatException e) {
		}
		return null;
	}

	private Double getXCorr(PSM psm) {
		try {
			if (psm != null) {
				final Set<edu.scripps.yates.utilities.proteomicsmodel.Score> scores = psm.getScores();
				for (final edu.scripps.yates.utilities.proteomicsmodel.Score score : scores) {
					if (score.getScoreName().toLowerCase().contains("xcorr")) {
						return Double.valueOf(score.getValue());
					}
				}
			}
		} catch (final NumberFormatException e) {
		}
		return null;
	}

	private CvParam getDeltaCnScore(Double deltacn) throws IOException {
		if (deltacn == null) {
			return null;
		}
		if (dtaSelectParser.getSearchEngines().contains(DTASelectParser.PROLUCID)) {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(DELTACN_PROLUCID,
					Score.getInstance(cvManager));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, String.valueOf(deltacn));
			return cvParam.getCvParam();
		} else {
			// if
			// (dtaSelectParser.getSearchEngines().contains(DTASelectParser.SEQUEST))
			// {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(DELTACN_SEQUEST,
					Score.getInstance(cvManager));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, String.valueOf(deltacn));
			return cvParam.getCvParam();
		}
		// return null;
	}

	private CvParam getXCorrScore(Double xcorr) throws IOException {
		if (xcorr == null) {
			return null;
		}
		final ControlVocabularyManager om = DTASelect2MzIdUtil.getOntologyManager();

		if (dtaSelectParser.getSearchEngines().contains(DTASelectParser.PROLUCID)) {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(XCORR_PROLUCID,
					Score.getInstance(om));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, String.valueOf(xcorr));
			return cvParam.getCvParam();
		} else {
			// if
			// (dtaSelectParser.getSearchEngines().contains(DTASelectParser.SEQUEST))
			// {
			final ControlVocabularyTerm cvTerm = DTASelect2MzIdUtil.getCvTermByAcc(XCORR_SEQUEST,
					Score.getInstance(om));
			final Param cvParam = DTASelect2MzIdUtil.getCVParam(cvTerm, String.valueOf(xcorr));
			return cvParam.getCvParam();
		}
		// return null;
	}

	private PeptideEvidenceRef getPeptideEvidenceRef(PSM dtaSelectPSM, Protein dtaSelectProtein) throws IOException {
		final PeptideEvidenceRef ret = new PeptideEvidenceRef();
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

			fileID = fileID.substring(1);
		}

		if (spectraDataBySpectraFileName.containsKey(fileID)) {
			return spectraDataBySpectraFileName.get(fileID);
		}

		final SpectraData spectraData = new SpectraData();
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
		final SpectrumIDFormat ret = new SpectrumIDFormat();
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
		final FileFormat format = new FileFormat();
		format.setCvParam(DTASelect2MzIdUtil
				.getCVParam(MS2_FORMAT_CV, "MS2 format", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		return format;
	}

	private FileFormat getMzXMLFileFormat() {
		final FileFormat format = new FileFormat();
		format.setCvParam(DTASelect2MzIdUtil
				.getCVParam(MzXML_FORMAT_CV, "ISB mzXML format", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		return format;
	}

	private FileFormat getDtaSelectFileFormat() {
		final FileFormat format = new FileFormat();
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
				final ParamList paramList = new ParamList();
				paramList.getCvParam().add(DTASelect2MzIdUtil
						.getCVParam("MS:1001494", "no threshold", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
				pdp.setThreshold(paramList);
			}
		}
		return pdp;
	}

	private ParamList getProteinThreshold() throws IOException {
		if (dtaSelectParser.getCommandLineParameter().getParametersMap().containsKey("--fp")) {
			final ParamList ret = new ParamList();
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
			final AnalysisSoftware searchEngine = getSearchEngine();
			String suffix = "";
			if (searchEngine != null) {
				suffix = "_" + searchEngine.getId();
			}
			pdl.setId("PDL" + suffix);
			final List<ProteinGroup> proteinGroups = getGroups();
			int numPAG = 0;
			for (final ProteinGroup proteinGroup : proteinGroups) {
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
		final SequenceCollection ret = new SequenceCollection();
		final Map<String, Protein> proteinMap = dtaSelectParser.getProteinMap();
		if (proteinMap != null) {
			for (final Protein protein : proteinMap.values()) {
				final DBSequence dbSequence = getDBSequence(protein);
				// check if already there
				if (!ret.getDBSequence().contains(dbSequence)) {
					ret.getDBSequence().add(dbSequence);
				}
				final List<PSM> psMs = protein.getPSMs();
				for (final PSM dtaSelectPSM : psMs) {
					final PeptideEvidence peptideEvidence = getPeptideEvidence(dtaSelectPSM, protein, ret);
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
		return ret;
	}

	private PeptideEvidence getPeptideEvidence(PSM dtaSelectPSM, Protein protein, SequenceCollection sequenceCollection)
			throws IOException {
		final String key = getPeptideEvidenceKey(dtaSelectPSM, protein);
		if (peptideEvidencesByKey.containsKey(key)) {
			return peptideEvidencesByKey.get(key);
		}
		final PeptideEvidence peptideEvidence = new PeptideEvidence();
		peptideEvidencesByKey.put(key, peptideEvidence);
		if (sequenceCollection != null) {
			sequenceCollection.getPeptideEvidence().add(peptideEvidence);
		}

		final DBSequence dbSequence = getDBSequence(protein);
		peptideEvidence.setDBSequence(dbSequence);
		final Peptide peptide = getPeptide(dtaSelectPSM, sequenceCollection);
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
			for (final PTM modification : dtaSelectPSM.getPTMs()) {
				for (final PTMSite site : modification.getPTMSites()) {
					peptide.getModification().add(getModification(modification, site));
				}
			}
		}
		// fix modifications. They are treated as different masses in the AA
		// in the search of PROLucid
		final List<Modification> fixedModifications = getFixedModifications(dtaSelectPSM);
		peptide.getModification().addAll(fixedModifications);
		final String peptideID = getPeptideKey(peptide);
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
		final StringBuilder sb = new StringBuilder();
		final String peptideSequence = peptide.getPeptideSequence();
		final List<Modification> modifications = peptide.getModification();
		Collections.sort(modifications, new Comparator<Modification>() {

			@Override
			public int compare(Modification o1, Modification o2) {
				final int location1 = o1.getLocation() != null ? o1.getLocation() : -1;
				final int location2 = o2.getLocation() != null ? o2.getLocation() : -1;
				return Integer.compare(location1, location2);
			}
		});
		// n term is location=0
		for (final Modification modification : modifications) {
			if (modification.getLocation() != null && modification.getLocation() == 0) {
				sb.append("(" + myFormatter4digits.format(modification.getMonoisotopicMassDelta()) + ")");
			}
		}
		for (int index = 0; index < peptideSequence.length(); index++) {
			final int position = index + 1;
			sb.append(peptideSequence.charAt(index));
			for (final Modification modification : modifications) {
				if (modification.getLocation() != null && modification.getLocation() == position) {
					sb.append("(" + myFormatter4digits.format(modification.getMonoisotopicMassDelta()) + ")");
				}
			}
		}
		// c term is location=peptide length+1
		for (final Modification modification : modifications) {
			if (modification.getLocation() != null && modification.getLocation() == peptideSequence.length() + 1) {
				sb.append("(" + myFormatter4digits.format(modification.getMonoisotopicMassDelta()) + ")");
			}
		}
		return sb.toString();
	}

	private List<PTM> getSortedByName(Collection<PTM> ptms) {
		final List<PTM> list = new ArrayList<PTM>();
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
		final List<Modification> ret = new ArrayList<Modification>();
		if (dtaSelectPSM != null) {
			final LabeledSearchType lst = getLabeledSearchTypeByFileName(dtaSelectPSM.getMSRun().getRunId());
			final SearchXmlFile searchParameters = SearchParametersManager.getInstance().getSearchParameters(lst);
			if (searchParameters != null) {

				// N term fix
				if (searchParameters.getNTermStaticMod() != null) {
					final double massDiff = Double.valueOf(searchParameters.getNTermStaticMod());
					if (Double.compare(massDiff, 0.0) != 0) {
						ret.add(getModification(null, 0, massDiff, null, true, false));
					}
				}

				// C term fix
				if (searchParameters.getCTermStaticMod() != null) {
					final double massDiff = Double.valueOf(searchParameters.getCTermStaticMod());
					if (Double.compare(massDiff, 0.0) != 0) {
						ret.add(getModification(null, dtaSelectPSM.getSequence().length(), massDiff, null, false,
								true));
					}
				}

				// modifications fix
				if (searchParameters.getStaticmods() != null) {

					// loop over the peptide sequence
					final String sequence = dtaSelectPSM.getSequence();
					for (int index = 0; index < sequence.length(); index++) {
						final char aa = sequence.charAt(index);
						for (final String diffAndResidue : searchParameters.getStaticmods()) {

							final double massDiff = Double.valueOf(diffAndResidue.split(" ")[0]);
							final String residues = diffAndResidue.split(" ")[1];
							for (int index2 = 0; index2 < residues.length(); index2++) {
								final char aa2 = residues.charAt(index2);
								if (aa == aa2) {
									boolean cTerminal = false;
									if (index + 1 == sequence.length()) {
										cTerminal = true;
									}
									ret.add(getModification(residues, index + 1, massDiff, null, false, cTerminal));
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
		final PeptideModificationUtil peptideModUtil = new PeptideModificationUtil(modification, site);
		final Modification ret = new Modification();
		ret.setLocation(peptideModUtil.getPosition());
		ret.setMonoisotopicMassDelta(peptideModUtil.getMonoDelta());
		ret.setAvgMassDelta(peptideModUtil.getAvgDelta());
		ret.getResidues().add(peptideModUtil.getResidues());
		if (peptideModUtil.getAccession() != null) {
			ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam(peptideModUtil.getAccession(), peptideModUtil.getName(),
					null, DTASelect2MzIdUtil.getPsimodCv()).getCvParam());
		} else {

			ret.getCvParam().add(DTASelect2MzIdUtil.getCVParam("MS:1001460", "unknown modification",
					peptideModUtil.getName(), DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());

		}

		return ret;

	}

	private DBSequence getDBSequence(Protein dtaSelectProtein) throws IOException {
		// look into the map
		if (dbSequences.containsKey(dtaSelectProtein.getAccession())) {
			return dbSequences.get(dtaSelectProtein.getAccession());
		}
		final DBSequence dbSequence = new DBSequence();
		dbSequences.put(dtaSelectProtein.getAccession(), dbSequence);
		final Accession acc = FastaParser.getACC(dtaSelectProtein.getAccession());
		dbSequence.setAccession(acc.getAccession());
		dbSequence.setId(getDBSequenceID(acc.getAccession()));
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
			final String fastaPath = dtaSelectParser.getFastaPath();
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
		final FileFormat ret = new FileFormat();
		ret.setCvParam(DTASelect2MzIdUtil
				.getCVParam("MS:1001348", "FASTA format", null, DTASelect2MzIdUtil.getPSIMsCv()).getCvParam());
		return ret;
	}

	private AnalysisSampleCollection getAnalysisSampleCollection(boolean createSampleIfNoExists) {
		final AnalysisSampleCollection ret = new AnalysisSampleCollection();
		final Sample sample2 = getSample(createSampleIfNoExists);
		if (sample2 == null) {
			return null;
		}
		ret.getSample().add(sample2);
		return ret;
	}

	private AuditCollection getAuditCollection() {
		final AuditCollection ret = new AuditCollection();
		ret.getOrganization().add(getYatesLabContact());
		ret.getPersonOrOrganization().add(getScrippsOrganization().getOrganization());
		return ret;
	}

	private Provider getDefaultProvider() {
		final Provider ret = new Provider();
		ret.setContactRole(getYatesLabContactRole(DTASelect2MzIdUtil.getResearcherRole()));
		ret.setId("PROVID");
		return ret;
	}

	private AnalysisSoftwareList getAnalysisSoftwareList() throws IOException {
		final AnalysisSoftwareList ret = new AnalysisSoftwareList();

		ret.getAnalysisSoftware().add(getSearchEngine());
		ret.getAnalysisSoftware().add(getDTASelectSoftware());
		return ret;
	}

	private AnalysisSoftware getDTASelectSoftware() throws IOException {
		final AnalysisSoftware ret = new AnalysisSoftware();
		ret.setId("DTASelect");
		ret.setName("DTASelect analysis software");
		ret.setVersion(dtaSelectParser.getSoftwareVersion());
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
		if (isProLuCIDSearch()) {
			return getProLuCIDAnalysisSoftware();
		}
		if (isSequestSearch()) {
			return getSequestAnalysisSoftware();
		}
		final AnalysisSoftware searchEngineSofware = getOtherSearchEngineSoftware();
		if (searchEngineSofware != null) {
			return searchEngineSofware;
		}
		return null;
	}

	private AnalysisSoftware getOtherSearchEngineSoftware() throws IOException {
		if (dtaSelectParser != null) {
			final Set<String> searchEngines = dtaSelectParser.getSearchEngines();
			for (final String searchEngine : searchEngines) {
				if (!searchEngine.equals(DTASelectParser.SEQUEST) && !searchEngine.equals(DTASelectParser.PROLUCID)) {
					final AnalysisSoftware ret = new AnalysisSoftware();
					ret.setId(searchEngine);
					ret.setName(searchEngine);
					ret.setVersion(dtaSelectParser.getSearchEngineVersion());

					final Param searchEngineUserParam = DTASelect2MzIdUtil.getUserParam(searchEngine, null);
					ret.setSoftwareName(searchEngineUserParam);
					return ret;
				}
			}
		}
		return null;
	}

	private AnalysisSoftware getSequestAnalysisSoftware() throws IOException {
		final AnalysisSoftware ret = new AnalysisSoftware();
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
		final AnalysisSoftware ret = new AnalysisSoftware();
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
		final ContactRole ret = new ContactRole();
		ret.setContact(getYatesLabContact());
		ret.setRole(role);
		return ret;
	}

	private Organization getYatesLabContact() {
		final Organization ret = new Organization();
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
		final ParentOrganization ret = new ParentOrganization();
		final Organization scrippsOrg = new Organization();
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
			for (final String searchEngine : searchEngines) {
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
		final CommandLineParser parser = new BasicParser();
		try {
			final CommandLine cmd = parser.parse(options, args);
			if (cmd.hasOption("i")) {
				final String iOptionValue = cmd.getOptionValue("i");
				inputFile = new File(iOptionValue);
				final File parentFile = inputFile.getParentFile();
				if (parentFile == null || !inputFile.exists()) {
					inputFile = new File(System.getProperty("user.dir") + File.separator + iOptionValue);
				}

				if (!inputFile.exists()) {
					errorInParameters("Input file '-i " + iOptionValue + "' doesn't exist or is not found");
				}
			} else {
				errorInParameters("Input file is missing");
			}

			if (cmd.hasOption("d")) {
				decoyRegexp = cmd.getOptionValue("d");
				try {
					Pattern.compile(decoyRegexp);
				} catch (final PatternSyntaxException e) {
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

				} catch (final Exception e) {
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
				final String v = cmd.getOptionValue("v");
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
				} catch (final Exception e) {
					errorInParameters("'" + sky + "' is not a valid boolean value");
				}
				log.info("Option 'sky' = '" + sky + "' detected. Compatibility with skyline is '" + skylineCompatible
						+ "'");
			} else {
				log.info("Option 'sky' not detected. Compatibility with skyline is by default '" + skylineCompatible
						+ "'");
			}

			final List<File> inputFiles = getInputFiles(inputFile, inputFileName, recursiveInputFileSearch);
			if (uniqueOutput) {
				final DTASelect2MzId conversor = new DTASelect2MzId(inputFiles, new File(System.getProperty("user.dir")
						+ File.separator + FilenameUtils.getBaseName(inputFileName) + ".mzid"), version);
				conversor.setDecoyRegexp(decoyRegexp);
				conversor.setIgnoreSpectra(ignoreSpectra);
				conversor.setSkylineCompatible(skylineCompatible);
				conversor.setReferenceToSpectra(referenceToSpectra);
				conversor.convert();
			} else {
				for (final File file : inputFiles) {
					inputFileName = FilenameUtils.getBaseName(file.getAbsolutePath());
					final File outputMzIdentMLFile = new File(
							file.getParentFile().getAbsolutePath() + File.separator + inputFileName + ".mzid");
					log.info("Using file name: '" + inputFileName + "' for output file(s).");
					final DTASelect2MzId conversor = new DTASelect2MzId(file, outputMzIdentMLFile, version);
					conversor.setDecoyRegexp(decoyRegexp);
					conversor.setIgnoreSpectra(ignoreSpectra);
					conversor.setSkylineCompatible(skylineCompatible);
					conversor.setReferenceToSpectra(referenceToSpectra);

					conversor.convert();
				}
			}
		} catch (final ParseException e) {
			errorInParameters(e.getMessage());
		} catch (final Exception e) {
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
		final List<File> ret = new ArrayList<File>();
		if (inputFile.isFile()) {
			ret.add(inputFile);
		} else {
			if (!recursiveInputFileSearch) {
				final File file = new File(inputFile.getAbsolutePath() + File.separator + inputFileName);
				ret.add(file);
			} else {
				final List<File> files = findFiles(inputFile, inputFileName);
				ret.addAll(files);
			}
		}
		final String plural = ret.size() > 1 ? "s" : "";
		log.info(ret.size() + " file" + plural + " found:");
		int i = 1;
		for (final File file : ret) {
			log.info(i++ + "- " + file.getAbsolutePath());
		}
		return ret;
	}

	private static List<File> findFiles(File inputFolder, String inputFileName) {
		final List<File> ret = new ArrayList<File>();
		// look in current folder
		final File file = new File(inputFolder.getAbsolutePath() + File.separator + inputFileName);
		if (file.exists()) {
			ret.add(file);
		}
		// look in subfolders
		final File[] listFiles = inputFolder.listFiles();
		for (final File file2 : listFiles) {
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
		final HelpFormatter formatter = new HelpFormatter();
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
