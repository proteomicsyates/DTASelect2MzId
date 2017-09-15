package edu.scripps.yates.dtaselect2mzid.skyline;

import java.awt.Color;
import java.awt.EventQueue;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.TextArea;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingWorker;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.border.TitledBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultStyledDocument;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;
import javax.swing.text.StyledDocument;

import org.apache.commons.io.FilenameUtils;

import com.sun.java.swing.plaf.windows.WindowsLookAndFeel;

import edu.scripps.yates.utilities.util.Pair;

public class MzXMLFixerGUI implements PropertyChangeListener {

	private JFrame frmFixMzxmlTo;
	private List<File> inputFiles = new ArrayList<File>();
	private Map<File, SwingWorker<Void, Void>> fixers = new HashMap<File, SwingWorker<Void, Void>>();
	private File outputFolder;
	private File currentFolder;
	private TextArea textArea;
	private JTextField outputFolderText;
	private JTextPane textPane;
	private JProgressBar progressBar;
	private StyledDocument document = new DefaultStyledDocument();
	private JLabel lblCancelProcess;
	private JButton btnCancel;
	private JButton selectButton;
	private JButton selectOutputButton;
	private JButton runButton;
	private boolean cancelling;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					MzXMLFixerGUI window = new MzXMLFixerGUI();
					window.frmFixMzxmlTo.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public MzXMLFixerGUI() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		try {
			UIManager.setLookAndFeel(new WindowsLookAndFeel());
		} catch (UnsupportedLookAndFeelException e) {

		}
		frmFixMzxmlTo = new JFrame();

		frmFixMzxmlTo.setTitle("Make mzXML compatible with Skyline library builder");
		frmFixMzxmlTo.setBounds(100, 100, 663, 606);
		frmFixMzxmlTo.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		GridBagLayout gridBagLayout = new GridBagLayout();
		gridBagLayout.columnWidths = new int[] { 200, 200 };
		gridBagLayout.rowHeights = new int[] { 30, 60, 30, 30, 30, 100, 30 };
		gridBagLayout.columnWeights = new double[] { 1.0, 0.0 };
		gridBagLayout.rowWeights = new double[] { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
		frmFixMzxmlTo.getContentPane().setLayout(gridBagLayout);

		JLabel lblSelectTheMzxml = new JLabel("Select the mzXML files you want to fix:");
		GridBagConstraints gbc_lblSelectTheMzxml = new GridBagConstraints();
		gbc_lblSelectTheMzxml.ipady = 10;
		gbc_lblSelectTheMzxml.anchor = GridBagConstraints.EAST;
		gbc_lblSelectTheMzxml.insets = new Insets(0, 0, 5, 5);
		gbc_lblSelectTheMzxml.gridx = 0;
		gbc_lblSelectTheMzxml.gridy = 0;
		frmFixMzxmlTo.getContentPane().add(lblSelectTheMzxml, gbc_lblSelectTheMzxml);

		selectButton = new JButton("Select");
		selectButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					selectInputFiles();
				} catch (URISyntaxException e1) {
					showError(e1);
				}
			}
		});
		GridBagConstraints gbc_selectButton = new GridBagConstraints();
		gbc_selectButton.anchor = GridBagConstraints.WEST;
		gbc_selectButton.insets = new Insets(0, 0, 5, 0);
		gbc_selectButton.gridx = 1;
		gbc_selectButton.gridy = 0;
		frmFixMzxmlTo.getContentPane().add(selectButton, gbc_selectButton);

		textArea = new TextArea();
		textArea.setRows(3);
		textArea.setEditable(false);
		GridBagConstraints gbc_textArea = new GridBagConstraints();
		gbc_textArea.anchor = GridBagConstraints.NORTH;
		gbc_textArea.fill = GridBagConstraints.HORIZONTAL;
		gbc_textArea.insets = new Insets(0, 0, 5, 0);
		gbc_textArea.gridwidth = 2;
		gbc_textArea.gridx = 0;
		gbc_textArea.gridy = 1;
		frmFixMzxmlTo.getContentPane().add(textArea, gbc_textArea);

		JLabel lblSelectTheOutput = new JLabel("Select the output folder:");
		GridBagConstraints gbc_lblSelectTheOutput = new GridBagConstraints();
		gbc_lblSelectTheOutput.anchor = GridBagConstraints.EAST;
		gbc_lblSelectTheOutput.insets = new Insets(0, 0, 5, 5);
		gbc_lblSelectTheOutput.gridx = 0;
		gbc_lblSelectTheOutput.gridy = 2;
		frmFixMzxmlTo.getContentPane().add(lblSelectTheOutput, gbc_lblSelectTheOutput);

		selectOutputButton = new JButton("Select");
		selectOutputButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					selectOutputFolder();
				} catch (URISyntaxException e1) {
					showError(e1);
				}
			}
		});
		GridBagConstraints gbc_selectOutputButton = new GridBagConstraints();
		gbc_selectOutputButton.anchor = GridBagConstraints.WEST;
		gbc_selectOutputButton.insets = new Insets(0, 0, 5, 0);
		gbc_selectOutputButton.gridx = 1;
		gbc_selectOutputButton.gridy = 2;
		frmFixMzxmlTo.getContentPane().add(selectOutputButton, gbc_selectOutputButton);

		outputFolderText = new JTextField();
		outputFolderText.setEditable(false);
		GridBagConstraints gbc_outputFolderText = new GridBagConstraints();
		gbc_outputFolderText.gridwidth = 2;
		gbc_outputFolderText.insets = new Insets(0, 0, 5, 0);
		gbc_outputFolderText.fill = GridBagConstraints.HORIZONTAL;
		gbc_outputFolderText.gridx = 0;
		gbc_outputFolderText.gridy = 3;
		frmFixMzxmlTo.getContentPane().add(outputFolderText, gbc_outputFolderText);
		outputFolderText.setColumns(10);

		JLabel lblClickToStart = new JLabel("Click to start:");
		GridBagConstraints gbc_lblClickToStart = new GridBagConstraints();
		gbc_lblClickToStart.anchor = GridBagConstraints.EAST;
		gbc_lblClickToStart.insets = new Insets(0, 0, 5, 5);
		gbc_lblClickToStart.gridx = 0;
		gbc_lblClickToStart.gridy = 4;
		frmFixMzxmlTo.getContentPane().add(lblClickToStart, gbc_lblClickToStart);

		runButton = new JButton("START");
		runButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				start();
			}
		});
		GridBagConstraints gbc_runButton = new GridBagConstraints();
		gbc_runButton.insets = new Insets(0, 0, 5, 0);
		gbc_runButton.anchor = GridBagConstraints.WEST;
		gbc_runButton.gridx = 1;
		gbc_runButton.gridy = 4;
		frmFixMzxmlTo.getContentPane().add(runButton, gbc_runButton);

		lblCancelProcess = new JLabel("Cancel process:");
		lblCancelProcess.setEnabled(false);
		lblCancelProcess.setToolTipText("Click here to cancel current process");
		GridBagConstraints gbc_lblCancelProcess = new GridBagConstraints();
		gbc_lblCancelProcess.anchor = GridBagConstraints.EAST;
		gbc_lblCancelProcess.insets = new Insets(0, 0, 5, 5);
		gbc_lblCancelProcess.gridx = 0;
		gbc_lblCancelProcess.gridy = 5;
		frmFixMzxmlTo.getContentPane().add(lblCancelProcess, gbc_lblCancelProcess);

		btnCancel = new JButton("Cancel");
		btnCancel.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				cancel();
			}
		});
		btnCancel.setEnabled(false);
		btnCancel.setToolTipText("Click here to cancel current process");
		GridBagConstraints gbc_btnCancel = new GridBagConstraints();
		gbc_btnCancel.anchor = GridBagConstraints.WEST;
		gbc_btnCancel.insets = new Insets(0, 0, 5, 0);
		gbc_btnCancel.gridx = 1;
		gbc_btnCancel.gridy = 5;
		frmFixMzxmlTo.getContentPane().add(btnCancel, gbc_btnCancel);

		textPane = new JTextPane(this.document);
		JScrollPane scroll = new JScrollPane(textPane);
		scroll.setViewportBorder(new TitledBorder(null, "Status", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		scroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		GridBagConstraints gbc_textPane = new GridBagConstraints();
		gbc_textPane.gridwidth = 2;
		gbc_textPane.insets = new Insets(0, 0, 5, 0);
		gbc_textPane.fill = GridBagConstraints.BOTH;
		gbc_textPane.gridx = 0;
		gbc_textPane.gridy = 6;
		frmFixMzxmlTo.getContentPane().add(scroll, gbc_textPane);

		progressBar = new JProgressBar();
		GridBagConstraints gbc_progressBar = new GridBagConstraints();
		gbc_progressBar.fill = GridBagConstraints.HORIZONTAL;
		gbc_progressBar.gridwidth = 2;
		gbc_progressBar.gridx = 0;
		gbc_progressBar.gridy = 7;
		frmFixMzxmlTo.getContentPane().add(progressBar, gbc_progressBar);

		showMessage("Welcome to the mzXML fixer!", true);
		showMessage(
				"With this program you will be able to fix your mzXML in order to use them into the Skyline library builder",
				false);
	}

	protected void cancel() {
		cancelling = true;

		for (File file : inputFiles) {
			if (fixers.containsKey(file)) {

				SwingWorker<Void, Void> swingWorker = fixers.get(file);
				showMessage("Trying to cancel processing of " + FilenameUtils.getName(file.getAbsolutePath()) + "...",
						true);

				if (!swingWorker.isDone() && !swingWorker.isCancelled()) {
					while (true) {
						boolean canceled = swingWorker.cancel(true);
						if (!canceled) {
							try {
								Thread.sleep(500);
							} catch (InterruptedException e) {
								e.printStackTrace();
							}
						} else {
							fixers.remove(file);
							break;
						}
					}
				}
			}
		}

		cancelling = false;
		enableControls(true);
	}

	protected void start() {
		checkInputs();
		if (inputFiles.isEmpty()) {
			if (textArea.getText().contains("\n")) {
				String[] split = textArea.getText().split("\n");
				for (String string : split) {
					inputFiles.add(new File(string));
				}
			}
		}
		for (File mzXMLFile : inputFiles) {
			try {

				MzXMLFixer fixer = new MzXMLFixer(mzXMLFile, outputFolder);
				fixer.addPropertyChangeListener(this);
				fixer.execute();
				fixers.put(mzXMLFile, fixer);
			} catch (IllegalArgumentException e) {
				showError(e);
			}
		}

	}

	private void checkInputs() {
		// TODO Auto-generated method stub

	}

	protected void selectOutputFolder() throws URISyntaxException {
		JFileChooser chooser = new JFileChooser(getCurrentFolder());
		chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		chooser.setMultiSelectionEnabled(false);
		chooser.setFileFilter(new FileNameExtensionFilter("mzXML", "*.*"));
		chooser.setApproveButtonText("Select");
		chooser.setApproveButtonToolTipText("Select the output folder where the fixed mzXML will be created");
		if (chooser.showOpenDialog(frmFixMzxmlTo) == JFileChooser.APPROVE_OPTION) {
			this.outputFolder = chooser.getSelectedFile();
			this.outputFolderText.setText(outputFolder.getAbsolutePath());
		}

	}

	private File getCurrentFolder() throws URISyntaxException {
		if (currentFolder == null) {
			currentFolder = new File(
					MzXMLFixerGUI.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath());
		}
		return currentFolder;
	}

	protected void showError(Exception e1) {
		showError(e1.getMessage());
	}

	protected void showError(String error) {
		StyleContext context = new StyleContext();
		Style style = context.getStyle(StyleContext.DEFAULT_STYLE);

		StyleConstants.setBold(style, true);
		StyleConstants.setForeground(style, Color.RED);
		try {
			this.document.insertString(this.document.getLength(), error + "\n", style);
		} catch (BadLocationException e) {
			showErrorDialog(e);
		}
	}

	private void showErrorDialog(BadLocationException e) {
		JOptionPane.showMessageDialog(frmFixMzxmlTo, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
	}

	protected void selectInputFiles() throws URISyntaxException {
		JFileChooser chooser = new JFileChooser(getCurrentFolder());
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		chooser.setMultiSelectionEnabled(true);
		chooser.setFileFilter(new FileNameExtensionFilter("mzXML files", "mzXML"));
		chooser.setApproveButtonText("Select");
		chooser.setApproveButtonToolTipText("Select one or more mzXML files");
		if (chooser.showOpenDialog(frmFixMzxmlTo) == JFileChooser.APPROVE_OPTION) {
			File[] selectedFiles = chooser.getSelectedFiles();
			for (File file : selectedFiles) {
				addInputFile(file);
			}
		}

	}

	private void addInputFile(File file) {
		this.inputFiles.add(file);
		setCurrentFolder(file.getParentFile());
		String previousText = textArea.getText();
		if (!"".equals(previousText)) {
			textArea.append("\n");
		}
		textArea.append(file.getAbsolutePath());

	}

	private void setCurrentFolder(File folder) {
		this.currentFolder = folder;

	}

	@Override
	public void propertyChange(PropertyChangeEvent evt) {
		File file = null;
		String message = null;
		switch (evt.getPropertyName()) {
		case MzXMLFixer.CANCELED:
			file = (File) evt.getNewValue();
			showError("Processing " + FilenameUtils.getName(file.getAbsolutePath()) + " cancelled");
			this.progressBar.setIndeterminate(false);

			break;
		case MzXMLFixer.DONE:
			file = (File) evt.getNewValue();
			while (cancelling) {
				try {
					Thread.sleep(500);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			inputFiles.remove(file);

			message = FilenameUtils.getName(file.getAbsolutePath()) + " fixed!";
			showMessage(message, true);
			if (inputFiles.isEmpty()) {
				enableControls(true);
				showMessage("All files fixed!", true);
				this.progressBar.setIndeterminate(false);
			}
			break;
		case MzXMLFixer.ERROR:
			this.progressBar.setIndeterminate(false);
			showError((Exception) evt.getNewValue());

			break;
		case MzXMLFixer.NUM_LINES:
			Pair<File, Long> pair = (Pair<File, Long>) evt.getNewValue();
			if (fixers.containsKey(pair.getFirstelement())) {
				Long numLines = pair.getSecondElement();
				if (numLines % 100 == 0) {
					showMessage(numLines + " lines fixed in "
							+ FilenameUtils.getName(pair.getFirstelement().getAbsolutePath()), false);
				}
			}
			break;
		case MzXMLFixer.STARTING:
			enableControls(false);

			file = (File) evt.getNewValue();
			message = "Processing " + FilenameUtils.getName(file.getAbsolutePath());
			showMessage(message, false);
			this.progressBar.setIndeterminate(true);
			break;
		default:
			break;
		}

	}

	private void showMessage(String message, boolean bold) {
		StyleContext context = new StyleContext();
		Style style = context.getStyle(StyleContext.DEFAULT_STYLE);
		if (bold) {
			StyleConstants.setBold(style, true);
		}
		try {
			this.document.insertString(this.document.getLength(), message + "\n", style);
		} catch (BadLocationException e) {
			showErrorDialog(e);
		}

	}

	private void enableControls(boolean enable) {
		this.btnCancel.setEnabled(!enable);
		this.lblCancelProcess.setEnabled(!enable);
		this.runButton.setEnabled(enable);
		this.selectButton.setEnabled(enable);
		this.selectOutputButton.setEnabled(enable);

	}
}
