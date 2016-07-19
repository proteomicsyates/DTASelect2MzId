package edu.scripps.yates.dtaselect2mzid.util;

public enum MzIdentMLVersion {
	VERSION_1_1("1.1.1"), VERSION_1_2("1.2.0");
	private final String versionString;

	private MzIdentMLVersion(String versionString) {
		this.versionString = versionString;
	}

	/**
	 * @return the versionString
	 */
	public String getVersionString() {
		return versionString;
	};

	public static MzIdentMLVersion getByString(String string) {
		for (MzIdentMLVersion version : values()) {
			if (version.getVersionString().equalsIgnoreCase(string)) {
				return version;
			}
		}
		return null;
	}

	public static String getCSVString() {
		StringBuilder sb = new StringBuilder();
		for (MzIdentMLVersion version : values()) {
			if (!"".equals(sb.toString())) {
				sb.append(", ");
			}
			sb.append("'" + version.getVersionString() + "'");
		}
		return sb.toString();
	}
}
