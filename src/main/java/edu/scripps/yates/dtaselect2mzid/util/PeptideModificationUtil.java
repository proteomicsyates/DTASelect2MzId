package edu.scripps.yates.dtaselect2mzid.util;

import java.net.URL;
import java.util.HashSet;
import java.util.Set;

import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import uk.ac.ebi.pridemod.PrideModController;
import uk.ac.ebi.pridemod.slimmod.model.SlimModCollection;
import uk.ac.ebi.pridemod.slimmod.model.SlimModification;

public class PeptideModificationUtil {
	private final PTM ptm;
	private final PTMSite ptmSite;
	private final SlimModification slimModification;
	private Double delta;
	private String residues;
	private static SlimModCollection preferredModifications;
	private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(PeptideModificationUtil.class);
	private static Set<String> errorMessages = new HashSet<String>();

	public PeptideModificationUtil(PTM dtaSelectModification, PTMSite ptmSite) {
		ptm = dtaSelectModification;
		this.ptmSite = ptmSite;
		if (preferredModifications == null) {
			URL url = getClass().getClassLoader().getResource("modification_mappings_dtaSelect.xml");
			if (url != null) {
				preferredModifications = PrideModController.parseSlimModCollection(url);
			} else {
				throw new IllegalArgumentException("Could not find preferred modification file");
			}
		}
		double delta = dtaSelectModification.getMassShift();
		double precision = 0.01;
		// map by delta
		SlimModCollection filteredMods = preferredModifications.getbyDelta(delta, precision);
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

	public PeptideModificationUtil(double massShift, String residues) {
		ptm = null;
		ptmSite = null;
		if (preferredModifications == null) {
			URL url = getClass().getClassLoader().getResource("modification_mappings_dtaSelect.xml");
			if (url != null) {
				preferredModifications = PrideModController.parseSlimModCollection(url);
			} else {
				throw new IllegalArgumentException("Could not find preferred modification file");
			}
		}
		this.residues = residues;
		delta = massShift;
		double precision = 0.01;
		// map by delta
		SlimModCollection filteredMods = preferredModifications.getbyDelta(delta, precision);
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

	public String getName() {
		if (slimModification != null)
			return slimModification.getShortNamePsiMod();
		return "unknown";
	}

	public String getAccession() {
		if (slimModification != null) {
			return slimModification.getIdPsiMod();
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
}
