package edu.scripps.yates.dtaselect2mzid.util;

import org.proteored.miapeapi.cv.ControlVocabularyManager;
import org.proteored.miapeapi.cv.ControlVocabularySet;

public class PeptideResultsDetailsCVSet extends ControlVocabularySet {

	private static PeptideResultsDetailsCVSet instance;

	public static PeptideResultsDetailsCVSet getInstance(ControlVocabularyManager cvManager) {
		if (instance == null)
			instance = new PeptideResultsDetailsCVSet(cvManager);
		return instance;
	}

	private PeptideResultsDetailsCVSet(ControlVocabularyManager cvManager) {
		super(cvManager);
		String[] parentAccessionsTMP = { "MS:1001105" };
		parentAccessions = parentAccessionsTMP;
		String[] explicitAccessionsTMP = {}; // slomo
												// score,
												// pepsplice
		explicitAccessions = explicitAccessionsTMP;

		miapeSection = 308;

	}

}
