package edu.scripps.yates.dtaselect2mzid.util;

public enum ReferenceToSpectra {
	MZXML, MS2;
	public static String getCSVString() {
		StringBuilder sb = new StringBuilder();
		for (ReferenceToSpectra referenceToSpectra : values()) {
			if (!"".equals(sb.toString())) {
				sb.append(", ");
			}
			sb.append("'" + referenceToSpectra.name() + "'");
		}
		return sb.toString();
	}
}
