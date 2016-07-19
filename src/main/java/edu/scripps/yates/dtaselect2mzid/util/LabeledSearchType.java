package edu.scripps.yates.dtaselect2mzid.util;

public enum LabeledSearchType {
	LIGHT(""), MEDIUM("M"), HEAVY("H");
	private final String key;

	private LabeledSearchType(String key) {
		this.key = key;
	}

	public String getKey() {
		return key;
	}
}
