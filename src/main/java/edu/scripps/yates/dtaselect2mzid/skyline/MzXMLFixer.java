package edu.scripps.yates.dtaselect2mzid.skyline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.SwingWorker;

import org.apache.commons.io.FilenameUtils;

import edu.scripps.yates.utilities.util.Pair;

/**
 * This class reads an input mzXML file and outputs the same file correcting
 * mz-int tag by m/z-int
 * 
 * @author Salva
 *
 */
public class MzXMLFixer extends SwingWorker<Void, Void> {
	private final File mzXMLFile;
	private final File outputFolder;
	public final static String STARTING = "MzXMLFixerSTARTING";
	public final static String CANCELED = "MzXMLFixerCANCELED";
	public final static String DONE = "MzXMLFixerDONE";
	public static final String ERROR = "MzXMLFixerERROR";
	public static final String NUM_LINES = "MzXMLFixerNUM_LINES";

	public MzXMLFixer(File mzXMLFile, File outputFolder) {
		if (mzXMLFile.getParentFile().getAbsolutePath().equals(outputFolder.getAbsolutePath())) {
			throw new IllegalArgumentException(
					"Output folder has to be different than the one containing the input mzXML file");
		}
		this.mzXMLFile = mzXMLFile;
		this.outputFolder = outputFolder;
	}

	public void fixMZXML() throws IOException {
		firePropertyChange(STARTING, null, this.mzXMLFile);
		File mzXML = this.mzXMLFile;
		File output = new File(
				this.outputFolder + File.separator + FilenameUtils.getName(this.mzXMLFile.getAbsolutePath()));
		FileWriter fw = null;
		BufferedReader br = null;
		try {
			long numLinesFixed = 0;
			fw = new FileWriter(output);
			br = new BufferedReader(new FileReader(mzXML));
			String line = br.readLine();
			while (line != null) {

				Thread.sleep(0);
				if (line.contains("\"mz-int\"")) {
					line = line.replace("\"mz-int\"", "\"m/z-int\"");
					numLinesFixed++;
					firePropertyChange(NUM_LINES, null, new Pair<File, Long>(this.mzXMLFile, numLinesFixed));
				}
				fw.write(line);
				fw.write("\n");
				line = br.readLine();
			}
		} catch (InterruptedException e) {

		} finally {
			try {
				fw.close();
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

		}
	}

	@Override
	protected Void doInBackground() throws Exception {
		try {
			fixMZXML();
		} catch (Exception e) {
			e.printStackTrace();
			firePropertyChange(ERROR, null, e);
		}
		return null;
	}

	@Override
	protected void done() {
		if (isCancelled()) {
			firePropertyChange(CANCELED, null, this.mzXMLFile);
		} else {
			firePropertyChange(DONE, null, this.mzXMLFile);
		}
		super.done();
	}
}
