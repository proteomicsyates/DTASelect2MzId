package edu.scripps.yates.dtaselect2mzid.util;

import java.net.URL;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import uk.ac.ebi.pride.utilities.pridemod.ModReader;
import uk.ac.ebi.pride.utilities.pridemod.model.Specificity;
import uk.ac.ebi.pride.utilities.pridemod.model.Specificity.AminoAcid;
import uk.ac.ebi.pridemod.PrideModController;
import uk.ac.ebi.pridemod.slimmod.model.SlimModCollection;
import uk.ac.ebi.pridemod.slimmod.model.SlimModification;

public class PeptideModificationUtil {
	private final PTM ptm;
	private final PTMSite ptmSite;
	private final SlimModification slimModification;
	private Double delta;
	private String residues;
	private final ModReader modReader;
	private uk.ac.ebi.pride.utilities.pridemod.model.PTM recognizedPTM;
	private static SlimModCollection preferredModifications;
	private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(PeptideModificationUtil.class);
	private static final String MOD0 = "MOD:00000";
	private static Set<String> errorMessages = new HashSet<String>();

	public PeptideModificationUtil(PTM modificaton, PTMSite ptmSite) {

		modReader = ModReader.getInstance();

		ptm = modificaton;
		this.ptmSite = ptmSite;
		if (preferredModifications == null) {
			final URL url = getClass().getClassLoader().getResource("modification_mappings_dtaSelect.xml");
			if (url != null) {
				preferredModifications = PrideModController.parseSlimModCollection(url);
			} else {
				throw new IllegalArgumentException("Could not find preferred modification file");
			}
		}
		final double delta = modificaton.getMassShift();
		final double precision = 0.01;
		// try with modreader
		final List<uk.ac.ebi.pride.utilities.pridemod.model.PTM> ptms = getPTMs(delta, residues);
		if (!ptms.isEmpty()) {
			recognizedPTM = ptms.get(0);
			residues = getResiduesString(recognizedPTM.getSpecificityCollection());
			slimModification = null;
		} else {
			// map by delta
			final SlimModCollection filteredMods = preferredModifications.getbyDelta(delta, precision);
			if (!filteredMods.isEmpty()) {
				slimModification = filteredMods.get(0);
			} else {
				final String message = "Peptide modification with delta mass=" + delta
						+ " is not recognized in the system. Please, contact system administrator in order to add it as a supported PTM in the system.";
				if (!errorMessages.contains(message)) {
					log.warn(message);
					errorMessages.add(message);
				}
				slimModification = null;
			}
		}
	}

	public PeptideModificationUtil(double massShift, String residues) {
		modReader = ModReader.getInstance();
		ptm = null;
		ptmSite = null;
		if (preferredModifications == null) {
			final URL url = getClass().getClassLoader().getResource("modification_mappings_dtaSelect.xml");
			if (url != null) {
				preferredModifications = PrideModController.parseSlimModCollection(url);
			} else {
				throw new IllegalArgumentException("Could not find preferred modification file");
			}
		}
		this.residues = residues;
		delta = massShift;

		// try with modreader
		final List<uk.ac.ebi.pride.utilities.pridemod.model.PTM> ptms = getPTMs(massShift, residues);
		if (!ptms.isEmpty()) {
			recognizedPTM = ptms.get(0);
			delta = recognizedPTM.getMonoDeltaMass();
			this.residues = getResiduesString(recognizedPTM.getSpecificityCollection());
			slimModification = null;
		} else {

			final double precision = 0.03;
			// map by delta
			final SlimModCollection filteredMods = preferredModifications.getbyDelta(delta, precision);
			if (!filteredMods.isEmpty()) {
				slimModification = filteredMods.get(0);
			} else {
				final String message = "Peptide modification with delta mass=" + delta
						+ " is not recognized in the system. Please, contact system administrator in order to add it as a supported PTM in the system.";
				if (!errorMessages.contains(message)) {
					log.warn(message);
					errorMessages.add(message);
				}
				slimModification = null;
			}
		}
	}

	private String getResiduesString(List<Specificity> specificityCollection) {
		final StringBuilder sb = new StringBuilder();
		for (final Specificity specificity : specificityCollection) {
			if (specificity.getName() == AminoAcid.NONE) {
				return "";
			}
			sb.append(specificity.getName());
		}
		return sb.toString();
	}

	private List<uk.ac.ebi.pride.utilities.pridemod.model.PTM> getPTMs(double massShift, String residues2) {
		final List<uk.ac.ebi.pride.utilities.pridemod.model.PTM> ptms = modReader.getPTMListByMonoDeltaMass(massShift);
		if (residues2 == null || "".equals(residues2)) {
			return ptms;
		}
		final List<uk.ac.ebi.pride.utilities.pridemod.model.PTM> ret = new ArrayList<uk.ac.ebi.pride.utilities.pridemod.model.PTM>();
		for (final uk.ac.ebi.pride.utilities.pridemod.model.PTM ptm : ptms) {
			final List<Specificity> specificityCollection = ptm.getSpecificityCollection();
			if (specificityCollection != null) {
				for (final Specificity specificity : specificityCollection) {
					final String aa = specificity.getName().toString();
					if (residues2.equals(aa)) {
						ret.add(ptm);
					}
				}
			}
		}
		if (ret.isEmpty()) {
			// if not found yet
			for (final uk.ac.ebi.pride.utilities.pridemod.model.PTM ptm : ptms) {
				final List<Specificity> specificityCollection = ptm.getSpecificityCollection();
				if (specificityCollection != null) {
					for (final Specificity specificity : specificityCollection) {
						if (specificity.getName() == AminoAcid.NONE) {
							ret.add(ptm);
						}
					}
				}
			}
		}
		return ret;
	}

	public String getName() {
		if (slimModification != null)
			return slimModification.getShortNamePsiMod();
		if (recognizedPTM != null) {
			return recognizedPTM.getName();
		}
		return "unknown";
	}

	public String getAccession() {
		if (slimModification != null) {
			final String idPsiMod = slimModification.getIdPsiMod();
			if (MOD0.equals(idPsiMod)) {
				return null;
			}
			return idPsiMod;
		}
		if (recognizedPTM != null) {
			return recognizedPTM.getAccession();
		}
		return null;
	}

	public int getPosition() {
		if (ptmSite != null) {
			return ptmSite.getPosition();
		}

		return -1;
	}

	public String getResidues() {
		if (ptmSite != null) {
			return ptmSite.getAA();
		}
		if (residues != null) {
			return residues;
		}
		return null;
	}

	public Double getMonoDelta() {
		if (ptm != null) {
			return ptm.getMassShift();
		}
		if (delta != null) {
			return delta;
		}
		return null;
	}

	public static Set<String> getErrorMessages() {
		return errorMessages;
	}

	public Double getAvgDelta() {
		if (recognizedPTM != null) {
			return recognizedPTM.getAveDeltaMass();
		}
		return null;
	}
}
