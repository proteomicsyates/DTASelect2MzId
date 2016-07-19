package edu.scripps.yates.dtaselect2mzid.util;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import uk.ac.ebi.pridemod.ModReader;

public class PeptideModificationUtil2 {
	private final PTM ptm;
	private final PTMSite ptmSite;
	private final ModReader modReader = ModReader.getInstance();
	private Double delta;
	private String residues;
	private static final org.apache.log4j.Logger log = org.apache.log4j.Logger
			.getLogger(PeptideModificationUtil2.class);
	private static Set<String> errorMessages = new HashSet<String>();
	private final uk.ac.ebi.pridemod.model.PTM pridePTM;

	public PeptideModificationUtil2(PTM dtaSelectModification, PTMSite ptmSite) {
		ptm = dtaSelectModification;
		this.ptmSite = ptmSite;

		double delta = dtaSelectModification.getMassShift();
		double precision = 0.01;
		// map by delta
		List<uk.ac.ebi.pridemod.model.PTM> pridePTMs = modReader.getPTMListByMonoDeltaMass(delta);
		if (pridePTMs != null && !pridePTMs.isEmpty()) {
			pridePTM = pridePTMs.get(0);
		} else {
			final String message = "Peptide modification with delta mass=" + delta
					+ " is not recognized in the system. Please, contact system administrator in order to add it as a supported PTM in the system.";
			if (!errorMessages.contains(message)) {
				log.warn(message);
				errorMessages.add(message);
			}
			pridePTM = null;
		}
	}

	public PeptideModificationUtil2(double massShift, String residues) {
		ptm = null;
		ptmSite = null;

		this.residues = residues;
		delta = massShift;
		double precision = 0.01;
		// map by delta
		List<uk.ac.ebi.pridemod.model.PTM> pridePTMs = modReader.getPTMListByMonoDeltaMass(delta);
		if (pridePTMs != null && !pridePTMs.isEmpty()) {
			pridePTM = pridePTMs.get(0);
		} else {
			final String message = "Peptide modification with delta mass=" + delta
					+ " is not recognized in the system. Please, contact system administrator in order to add it as a supported PTM in the system.";
			if (!errorMessages.contains(message)) {
				log.warn(message);
				errorMessages.add(message);
			}
			pridePTM = null;
		}
	}

	public String getName() {
		if (pridePTM != null)
			return pridePTM.getName();
		return "unknown";
	}

	public String getAccession() {
		if (pridePTM != null) {
			return pridePTM.getAccession();
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
