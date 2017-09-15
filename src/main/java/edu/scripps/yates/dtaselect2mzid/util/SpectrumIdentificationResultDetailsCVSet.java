package edu.scripps.yates.dtaselect2mzid.util;

import org.proteored.miapeapi.cv.ControlVocabularyManager;
import org.proteored.miapeapi.cv.ControlVocabularySet;

public class SpectrumIdentificationResultDetailsCVSet extends ControlVocabularySet {

	private static SpectrumIdentificationResultDetailsCVSet instance;

	public static SpectrumIdentificationResultDetailsCVSet getInstance(ControlVocabularyManager cvManager) {
		if (instance == null)
			instance = new SpectrumIdentificationResultDetailsCVSet(cvManager);
		return instance;
	}

	private SpectrumIdentificationResultDetailsCVSet(ControlVocabularyManager cvManager) {
		super(cvManager);
		String[] parentAccessionsTMP = { "MS:1001405" };
		parentAccessions = parentAccessionsTMP;
		String[] explicitAccessionsTMP = {}; // slomo
												// score,
												// pepsplice
		explicitAccessions = explicitAccessionsTMP;

		miapeSection = 308;

	}

}
