package phylip;

import java.awt.EventQueue;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.Math;
import java.util.Scanner;

import javax.swing.JFrame;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.JLabel;
import javax.swing.SwingConstants;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JFileChooser;

import utilities.DrawPreview;

import java.awt.Font;
import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class DrawtreeUserInterface {


	public class DrawtreeData{
		String  intree;
		String  plotfile;
		String  plotfileopt;
		String  usefont;
		String  treegrows;
		boolean usebranchlengths;
		String  labeldirec;
		Double  labelangle;
		Double  treerotation;
		Double  treearc;
		String  iterationkind;
		int     iterationcount;
		boolean regularizeangles;
		boolean avoidlabeloverlap;
		boolean branchrescale;
		Double  branchscaler;
		Double  relcharhgt;
		Double  xmarginratio;
		Double  ymarginratio;
		String  librarypath;
		boolean doplot; 
		String  finalplottype;
	}

	private String filedir;
	private boolean phylipCall;
	private boolean plotCall;
	private String calledBy;
	
	private JFrame frmDrawtreeControls;
	private JButton btnInputTree;
	private JTextField txtInputTree;
	private JTextField txtPlot;
	private JButton btnPlot;
	private JComboBox cmbxPlotFont;
	private JRadioButton rbtnUseLenY;
	private JRadioButton rbtnUseLenN;
	private JComboBox cmbxLabelAngle;
	private JTextField txtTreeRotation;
	private JLabel lblTreeRotation;
	private JLabel lblFixedAngleOf;
	private JComboBox cmbxAngle;
	private JComboBox cmbxIterate;
	private JLabel lblAvoidOverlap;
	private JRadioButton rbtnAvoidOverY;
	private JRadioButton rbtnAvoidOverN;
	private JComboBox cmbxRescale;
	private JTextField txtRelCharHgt;
	private JLabel lblRegularizeTheAngles;
	private JRadioButton rbtnRegAngleY;
	private JRadioButton rbtnRegAngleN;
	private JLabel lblBranchScale;
	private JTextField txtBranchScale;
	private JTextField txtTreeArc;
	private JRadioButton rbtnTreeH;
	private JRadioButton rbtnTreeV;
	private JLabel lblMaximumIterations;
	private JTextField txtIterationCount;
	private JLabel lblFinalPlotType;
	private JComboBox cmbxFinalPlotType;
	private JLabel lblTreeArc;
	private JLabel lblIfPresent;
	private JButton btnPreview; 
	private JButton btnQuit; 
	private JLabel lblAngleOfLabels;
	private JLabel lblIterateToImprove;
	private JLabel lblBranchLengths;
	private JLabel lblGrows;
	private JLabel lblUseBranchLengths;
	private JLabel lblFont;
	private JButton btnPlotFile;
	private JLabel lblRelCharHgt;
	private JLabel lblMarginRatios;
	private JLabel lblXMarginRatio;
	private JTextField txtXMarginRatio;
	private JLabel lblYMarginRatio;
	private JTextField txtYMarginRatio;

	private JScrollPane scrollPane;
	private JPanel panel;
	private JButton btnRestoreDefaults;
	private JButton btnReadInit;

	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DrawtreeUserInterface window = new DrawtreeUserInterface(args);
					window.frmDrawtreeControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	protected void ChooseFile(JTextField file) {		
		//Construct a new file choose whose default path is the path to this executable, which 
		//is returned by System.getProperty("user.dir")
		
		JFileChooser fileChooser = new JFileChooser( filedir);

		int option = fileChooser.showOpenDialog(frmDrawtreeControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}	
	}
	
	protected void BranchLengthToggle(boolean uselength) {	
		
		if (uselength){
			rbtnUseLenY.setSelected(true);
			rbtnUseLenN.setSelected(false);
		}	
		else{
			rbtnUseLenY.setSelected(false);
			rbtnUseLenN.setSelected(true);
		}		
	}
	
	protected void AvoidOverlapToggle(boolean avoidoverlap) {		
		if (avoidoverlap){
			rbtnAvoidOverY.setSelected(true);
			rbtnAvoidOverN.setSelected(false);
		}	
		else{
			rbtnAvoidOverY.setSelected(false);
			rbtnAvoidOverN.setSelected(true);
		}		
	}
	
	protected void LabelAngleToggle() {		
		if ((cmbxLabelAngle.getSelectedItem().toString()).contains("Fixed")){
			lblFixedAngleOf.setEnabled(true);
			cmbxAngle.setEnabled(true);
		}	
		else{
			lblFixedAngleOf.setEnabled(false);
			cmbxAngle.setEnabled(false);
		}		
	}
	
	protected void IterationType() {		
		if ((cmbxIterate.getSelectedItem().toString()).contains("No")){
			lblRegularizeTheAngles.setEnabled(true);
			rbtnRegAngleY.setEnabled(true);
			rbtnRegAngleN.setEnabled(true);
			lblAvoidOverlap.setEnabled(false);
			rbtnAvoidOverY.setEnabled(false);
			rbtnAvoidOverN.setEnabled(false);
			lblMaximumIterations.setEnabled(false);
			txtIterationCount.setEnabled(false);
			txtIterationCount.setText("0");
		}	
		else{
			lblRegularizeTheAngles.setEnabled(false);
			rbtnRegAngleY.setEnabled(false);
			rbtnRegAngleN.setEnabled(false);
			lblAvoidOverlap.setEnabled(true);
			rbtnAvoidOverY.setEnabled(true);
			rbtnAvoidOverN.setEnabled(true);
			lblMaximumIterations.setEnabled(true);
			txtIterationCount.setEnabled(true);
			if ((cmbxIterate.getSelectedItem().toString()).contains("Equal")){
				txtIterationCount.setText("100");
			}
			else{
				txtIterationCount.setText("50");
			}
		}		
	}
	
	protected void ScaleValue() {		
		if ((cmbxRescale.getSelectedItem().toString()).contains("Fixed")){
			lblBranchScale.setEnabled(true);
			txtBranchScale.setEnabled(true);
			txtBranchScale.setEditable(true);
		}	
		else{
			lblBranchScale.setEnabled(false);
			txtBranchScale.setEnabled(false);
			txtBranchScale.setEditable(false);
		}		
	}
	
	protected void RegAngleToggle(boolean doregular) {		
		if (doregular){
			rbtnRegAngleY.setSelected(true);
			rbtnRegAngleN.setSelected(false);
		}	
		else{
			rbtnRegAngleY.setSelected(false);
			rbtnRegAngleN.setSelected(true);
		}		
	}
	
	protected void SweepLimit() {
		// doing a mod 360 in case the user gets clever
		if (Math.abs(Double.parseDouble(txtTreeArc.getText())) > 360)
		{
			txtTreeArc.setText(Double.toString(Double.parseDouble(txtTreeArc.getText())%360));
		}
		if (Double.parseDouble(txtTreeArc.getText()) == 0.0)
		{
			txtTreeArc.setText(Double.toString(360));
		}
	}
	
	protected void RotationLimit() {
		// doing a mod 360 in case the user gets clever
		if (Double.parseDouble(txtTreeRotation.getText()) > 360)
		{
			txtTreeRotation.setText(Double.toString(Double.parseDouble(txtTreeRotation.getText())%360));
		}
	}
	
	protected void TreeGrowToggle(boolean ishoriz) {		
		if (ishoriz){
			rbtnTreeH.setSelected(true);
			rbtnTreeV.setSelected(false);
		}	
		else{
			rbtnTreeH.setSelected(false);
			rbtnTreeV.setSelected(true);
		}		
	}

	protected boolean LaunchDrawtreeInterface(Boolean doPlot){
		DrawtreeData inputdata = getInputVals();
		inputdata.doplot = doPlot;	
		DrawtreeInterface dg = new DrawtreeInterface();
		return (dg.DrawtreeRun(inputdata));		
	}	

	/**
	 * Create the application.
	 */
	public DrawtreeUserInterface(String[] args) {
		phylipCall = false;
		plotCall = false;
		if (args.length > 0)
		{
			// currently the effect of both these calls is the same, 
			// but they could easily diverge in the future
			calledBy = args[0];
			if (args[0].contains("Phylip"))
			{
				// called from the Phylip super program
				phylipCall = true;
			}
			else
			{
				// called from some other Phylip program
				plotCall = true;				
			}
		}
		initialize();
		DoInit();
		
		if (plotCall){
			String sourceTree = calledBy + "Outtree";
  			txtInputTree.setText(sourceTree);
 		}
	}
	
	protected void DoInit()
	{
		// reset everything if there is an init file
		getStoredSettings();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		filedir = System.getProperty("user.dir");
		
		frmDrawtreeControls = new JFrame();
		frmDrawtreeControls.getContentPane().setBackground(new Color(204, 255, 255));
		frmDrawtreeControls.setBackground(new Color(204, 255, 255));
		frmDrawtreeControls.setTitle("Drawtree");
		frmDrawtreeControls.setBounds(100, 100, 650, 650);
		frmDrawtreeControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDrawtreeControls.setPreferredSize(new Dimension(frmDrawtreeControls.getBounds().width, frmDrawtreeControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmDrawtreeControls.getPreferredSize());
		frmDrawtreeControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmDrawtreeControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow]", "[][][][]"));
		
		btnInputTree = new JButton("Input Tree");
		btnInputTree.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputTree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtInputTree);
			}
		});
		panel.add(btnInputTree, "cell 0 0,growx");
		
		txtInputTree = new JTextField();
		txtInputTree.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInputTree, "cell 1 0 2 1,growx");
		
		btnPlot = new JButton("Plot File");
		btnPlot.setFont(new Font("Arial", Font.BOLD, 13));
		btnPlot.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtPlot);
			}
		});
		panel.add(btnPlot, "cell 0 1,growx");
		
		txtPlot = new JTextField();
		txtPlot.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtPlot, "cell 1 1 2 1,growx");
		
		lblFont = new JLabel("PostScript Font:");
		lblFont.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblFont, "flowx,cell 1 2 1 1,alignx right");
		
		cmbxPlotFont = new JComboBox();
		cmbxPlotFont.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxPlotFont.setModel(new DefaultComboBoxModel(new String[] {
				"Times-Roman", "Times-Bold", "Helvetica", "Helvetica-Bold", "Courier", "Courier-Bold", 
				"AvantGarde-Book","AvantGarde-BookOblique", "AvantGarde-Demi", "AvantGarde-DemiOblique", 
				"Bookman-Light", "Bookman-LightItalic", "Bookman-Demi", "Bookman-DemiItalic", 
				"Courier-Oblique", "Courier-BoldOblique", "Helvetica-Oblique", 
				"Helvetica-BoldOblique", "Helvetica-Narrow", "Helvetica-Narrow-Oblique", 
				"Helvetica-Narrow-Bold", "Helvetica-Narrow-BoldOblique", "NewCenturySchlbk-Roman", 
				"NewCenturySchlbk-Italic", "NewCenturySchlbk-Bold", "NewCenturySchlbk-BoldItalic", 
				"Palatino-Roman", "Palatino-Italic", "Palatino-Bold", "Palatino-BoldItalic", 
				 "Times-BoldItalic",  "Times-Italic", "ZapfChancery-MediumItalic"}));
		panel.add(cmbxPlotFont, "cell 2 2,growx");	
		
		lblGrows = new JLabel("Tree grows:");
		lblGrows.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGrows, "flowx,cell 1 3 1 1,alignx right");
		
		rbtnTreeH = new JRadioButton("Horizontally");
		rbtnTreeH.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnTreeH.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				TreeGrowToggle(true);
			}
		});
		panel.add(rbtnTreeH, "cell 2 3");	
		
		rbtnTreeV = new JRadioButton("Vertically");
		rbtnTreeV.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnTreeV.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				TreeGrowToggle(false);
			}
		});
		panel.add(rbtnTreeV, "cell 2 3");	
		
		lblUseBranchLengths = new JLabel("Use branch lengths:");
		lblUseBranchLengths.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseBranchLengths, "flowx,cell 1 4 1 1,alignx right");
		
		rbtnUseLenY = new JRadioButton("Yes");
		rbtnUseLenY.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnUseLenY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				BranchLengthToggle(true);
			}
		});
		panel.add(rbtnUseLenY, "cell 2 4");	
		
		rbtnUseLenN = new JRadioButton("No");
		rbtnUseLenN.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnUseLenN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				BranchLengthToggle(false);
			}
		});
		panel.add(rbtnUseLenN, "cell 2 4");	
		
		lblIfPresent = new JLabel(" (if present)");
		lblIfPresent.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblIfPresent, "cell 2 4");	
		
		lblAngleOfLabels = new JLabel("Angle of labels:");
		lblAngleOfLabels.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAngleOfLabels, "flowx,cell 1 5 1 1,alignx right");
		
		cmbxLabelAngle = new JComboBox();
		cmbxLabelAngle.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxLabelAngle.setModel(new DefaultComboBoxModel(new String[] {"Middle of Label", "Fixed", "Radial", "Along Branches"}));
		cmbxLabelAngle.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				LabelAngleToggle();
			}
		});
		panel.add(cmbxLabelAngle, "cell 2 5,growx");	
		
		lblFixedAngleOf = new JLabel("Fixed label angle:");
		lblFixedAngleOf.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblFixedAngleOf, "flowx,cell 1 6 1 1,alignx right");
		
		cmbxAngle = new JComboBox();
		cmbxAngle.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxAngle.setModel(new DefaultComboBoxModel(new String[] {"    0.0", "  90.0", "-90.0"}));
		panel.add(cmbxAngle, "cell 2 6");	

		lblTreeRotation = new JLabel("Angle of tree:");
		lblTreeRotation.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTreeRotation, "flowx,cell 1 7 1 1,alignx right");
		
		txtTreeRotation = new JTextField();
		txtTreeRotation.setFont(new Font("Arial", Font.PLAIN, 13));
		txtTreeRotation.setHorizontalAlignment(SwingConstants.RIGHT);
		txtTreeRotation.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				RotationLimit();
			}
		});
		txtTreeRotation.setColumns(6);
		panel.add(txtTreeRotation, "cell 2 7");	
		
		lblTreeArc = new JLabel("Arc of tree:");
		lblTreeArc.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTreeArc, "flowx,cell 1 8 1 1,alignx right");
		
		txtTreeArc = new JTextField();
		txtTreeArc.setFont(new Font("Arial", Font.PLAIN, 13));
		txtTreeArc.setHorizontalAlignment(SwingConstants.RIGHT);
		txtTreeArc.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				SweepLimit();
			}
		});
		txtTreeArc.setColumns(6);
		panel.add(txtTreeArc, "cell 2 8");	
		
		lblIterateToImprove = new JLabel("Iterate to improve tree:");
		lblIterateToImprove.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblIterateToImprove, "flowx,cell 0 9 2 1,alignx right");
		
		cmbxIterate = new JComboBox();
		cmbxIterate.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxIterate.setModel(new DefaultComboBoxModel(new String[] {"Equal-Daylight algorithm", "n-Body algorithm", "No"}));
		cmbxIterate.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				IterationType();
			}
		});
		panel.add(cmbxIterate, "cell 2 9,growx");	
		
		lblMaximumIterations = new JLabel("Maximum iterations:");
		lblMaximumIterations.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMaximumIterations, "flowx,cell 1 10 1 1,alignx right");
		
		txtIterationCount = new JTextField();
		txtIterationCount.setFont(new Font("Arial", Font.PLAIN, 13));
		txtIterationCount.setColumns(6);
		panel.add(txtIterationCount, "cell 2 10");	
		
		lblRegularizeTheAngles = new JLabel("Regularize the angles:");
		lblRegularizeTheAngles.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRegularizeTheAngles, "flowx,cell 1 11 1 1,alignx right");
		
		rbtnRegAngleY = new JRadioButton("Yes");
		rbtnRegAngleY.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRegAngleY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				RegAngleToggle(true);
			}
		});
		panel.add(rbtnRegAngleY, "cell 2 11");	
		
		rbtnRegAngleN = new JRadioButton("No");
		rbtnRegAngleN.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRegAngleN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				RegAngleToggle(false);
			}
		});
		panel.add(rbtnRegAngleN, "cell 2 11");	
		
		lblAvoidOverlap = new JLabel("Try to avoid label overlap:");
		lblAvoidOverlap.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAvoidOverlap, "flowx,cell 1 12 1 1,alignx right");
		
		rbtnAvoidOverY = new JRadioButton("Yes");
		rbtnAvoidOverY.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAvoidOverY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				AvoidOverlapToggle(true);
			}
		});
		panel.add(rbtnAvoidOverY, "cell 2 12");	
		
		rbtnAvoidOverN = new JRadioButton("No");
		rbtnAvoidOverN.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAvoidOverN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				AvoidOverlapToggle(false);
			}
		});
		panel.add(rbtnAvoidOverN, "cell 2 12");	
		
		lblBranchLengths = new JLabel("Branch lengths:");
		lblBranchLengths.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBranchLengths, "flowx,cell 1 13 1 1,alignx right");
		
		cmbxRescale = new JComboBox();
		cmbxRescale.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxRescale.setModel(new DefaultComboBoxModel(new String[] {"Automatically rescale", "Fixed scale"}));
		cmbxRescale.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ScaleValue();
			}
		});
		panel.add(cmbxRescale, "cell 2 13,growx");	
		
		lblRelCharHgt = new JLabel("Relative character height:");
		lblRelCharHgt.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRelCharHgt, "flowx,cell 1 15 1 1,alignx right");
		
		txtRelCharHgt = new JTextField();
		txtRelCharHgt.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRelCharHgt.setColumns(6);
		panel.add(txtRelCharHgt, "cell 2 15");	
				
		lblBranchScale = new JLabel("Branch scale:");
		lblBranchScale.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBranchScale, "flowx,cell 1 14 1 1,alignx right");
		
		txtBranchScale = new JTextField();
		txtBranchScale.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBranchScale.setEditable(false);
		txtBranchScale.setColumns(6);
		panel.add(txtBranchScale, "cell 2 14");	
		
		lblFinalPlotType = new JLabel("Final plot file type:");
		lblFinalPlotType.setHorizontalAlignment(SwingConstants.TRAILING);
		lblFinalPlotType.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblFinalPlotType, "flowx,cell 1 16 1 1,alignx right");
		
		cmbxFinalPlotType = new JComboBox();
		cmbxFinalPlotType.setModel(new DefaultComboBoxModel(new String[] {"Postscript", "SVG"}));
		cmbxFinalPlotType.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxFinalPlotType, "cell 2 16,growx");	
		
		lblMarginRatios = new JLabel("Page margin ratios:");
		lblMarginRatios.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMarginRatios.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMarginRatios, "flowx,cell 1 17 1 1,alignx right");
		
		lblXMarginRatio = new JLabel("X:");
		lblXMarginRatio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblXMarginRatio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblXMarginRatio, "cell 2 17");	
		
		txtXMarginRatio = new JTextField();
		txtXMarginRatio.setColumns(6);
		panel.add(txtXMarginRatio, "cell 2 17");	
		
		lblYMarginRatio = new JLabel("Y:");
		lblYMarginRatio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblYMarginRatio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblYMarginRatio, "cell 2 17");	
		
		txtYMarginRatio = new JTextField();
		txtYMarginRatio.setColumns(6);
		panel.add(txtYMarginRatio, "cell 2 17");	
		
		btnRestoreDefaults = new JButton("Restore Defaults");
		btnRestoreDefaults.setFont(new Font("Arial", Font.BOLD, 13));
		btnRestoreDefaults.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				resetDefaults();
			}		
		});
		panel.add(btnRestoreDefaults, "cell 0 3,growx");
		
		btnReadInit = new JButton("Read Init file");
		btnReadInit.setFont(new Font("Arial", Font.BOLD, 13));
		btnReadInit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				getStoredSettings();
			}		
		});
		panel.add(btnReadInit, "cell 0 2,growx");
	
		btnPreview = new JButton("Preview");
		btnPreview.setFont(new Font("Arial", Font.BOLD, 13));
		btnPreview.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean retval = LaunchDrawtreeInterface(false);
				if (retval)
				{
					String title = "Preview: " + (String)txtPlot.getText();
					String curDir = System.getProperty("user.dir");
					curDir += "/JavaPreview.ps";
					new DrawPreview(title, curDir);
				}
			}
		});
		panel.add(btnPreview, "cell 2 18,alignx right");	
		
		btnPlotFile = new JButton("Plot");
		btnPlotFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnPlotFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				@SuppressWarnings("unused")
				boolean retval = LaunchDrawtreeInterface(true);
			}
		});
		panel.add(btnPlotFile, "cell 2 18");	
		
		btnQuit = new JButton("Quit");
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveSettings();
				if((phylipCall) || (plotCall))
				{
					frmDrawtreeControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		panel.add(btnQuit, "cell 2 18");	
		
	}
	
	protected DrawtreeData getInputVals()
	{
		DrawtreeData inputdata = new DrawtreeData();
				
		inputdata.intree = (String)txtInputTree.getText();
		inputdata.plotfile = (String)txtPlot.getText();
		inputdata.usefont =  cmbxPlotFont.getSelectedItem().toString();	
		inputdata.usebranchlengths = rbtnUseLenY.isSelected();
		
		switch (cmbxLabelAngle.getSelectedIndex()){
			case 0: //middle
				inputdata.labeldirec = "middle";
				break;
			case 1: //fixed
				inputdata.labeldirec = "fixed";
				break;
			case 2: //radial
				inputdata.labeldirec = "radial";
				break;
			case 3: //along
				inputdata.labeldirec = "along";
				break;
			default:
				inputdata.labeldirec = "middle";
				break;													
		}
		
		switch (cmbxAngle.getSelectedIndex()){
			case 0: //0.0
				inputdata.labelangle = 0.0;
				break;
			case 1: //90.0
				inputdata.labelangle = 90.0;
				break;
			case 2: //-90.0
				inputdata.labelangle = -90.0;
				break;
			default:
				inputdata.labelangle = 0.0;
				break;													
		}
		
		if (rbtnTreeH.isSelected())
		{
			inputdata.treegrows = "horizontal";
		}
		else
		{
			inputdata.treegrows = "vertical";					
		}
		
		// just in case the user managed to enter values without a carriage return
		RotationLimit();
		SweepLimit();
		
		inputdata.treerotation = new Double(txtTreeRotation.getText());
		inputdata.treearc = new Double(txtTreeArc.getText());
		
		switch (cmbxIterate.getSelectedIndex()){
			case 0: //equal daylight
				inputdata.iterationkind = "improve";
				break;
			case 1: //nbody
				inputdata.iterationkind = "nbody";
				break;
			case 2: //no
				inputdata.iterationkind = "no";
				break;
			default:
				inputdata.iterationkind = "improve";
				break;													
		}
		
		inputdata.iterationcount = Integer.parseInt(txtIterationCount.getText());

		inputdata.regularizeangles = rbtnRegAngleY.isSelected();
		inputdata.avoidlabeloverlap = rbtnAvoidOverY.isSelected();
		
		switch (cmbxRescale.getSelectedIndex()){
			case 0: //automatic
				inputdata.branchrescale = true;
				break;
			case 1: //fixed length
				inputdata.branchrescale = false;
				break;
			default:
				inputdata.branchrescale = true;
				break;													
		}
		inputdata.branchscaler = new Double(txtBranchScale.getText());
		inputdata.relcharhgt = new Double(txtRelCharHgt.getText());
		inputdata.xmarginratio = new Double(txtXMarginRatio.getText());
		inputdata.ymarginratio = new Double(txtYMarginRatio.getText());
		inputdata.librarypath = System.getProperty("user.dir");
		
		switch (cmbxFinalPlotType.getSelectedIndex()){
			case 0: // postscript
				inputdata.finalplottype = "lw";
				break;
			case 1: // svg
				inputdata.finalplottype = "svg";
				break;
			default:
				inputdata.finalplottype = "lw";
				break;	
		}
		
		return inputdata;
		
	}

	
	protected void saveSettings(){
		DrawtreeData inputvals = getInputVals();
		// there must be a better way to format this output, but this works for the prototype JRM
        try {
            BufferedWriter output = new BufferedWriter(new FileWriter("drawtreeInit.txt"));
    		output.write("intree : "+inputvals.intree+"\n");
    		output.write("plotfile : "+inputvals.plotfile+"\n");
    		output.write("usefont : "+cmbxPlotFont.getSelectedIndex()+" : "+inputvals.usefont+"\n");
    		output.write("usebranchlengths : "+String.format("%b",inputvals.usebranchlengths)+"\n");		
    		output.write("labeldirec : "+cmbxLabelAngle.getSelectedIndex()+" : "+inputvals.labeldirec+"\n");	
            output.write("labelangle : "+cmbxAngle.getSelectedIndex()+" : "+String.format("%.1f",inputvals.labelangle)+"\n");
            output.write("treegrows : "+inputvals.treegrows+"\n");
    		output.write("treerotation : "+String.format("%.1f",inputvals.treerotation)+"\n");
    		output.write("treearc : "+String.format("%.1f",inputvals.treearc)+"\n");	
    		output.write("iterationkind : "+cmbxIterate.getSelectedIndex()+" : "+inputvals.iterationkind+"\n");	
    		output.write("iterationcount : "+inputvals.iterationcount+"\n");
    		output.write("regularizeangles : "+String.format("%b",inputvals.regularizeangles)+"\n");
    		output.write("avoidlabeloverlap : "+String.format("%b",inputvals.avoidlabeloverlap)+"\n");		
    		output.write("branchrescale : "+cmbxRescale.getSelectedIndex()+" : "+String.format("%b",inputvals.branchrescale)+"\n");
    		output.write("branchscaler : "+String.format("%.1f",inputvals.branchscaler)+"\n");
    		output.write("relcharhgt : "+String.format("%.4f",inputvals.relcharhgt)+"\n");
    		output.write("xmarginratio : "+String.format("%.1f",inputvals.xmarginratio)+"\n");
    		output.write("ymarginratio : "+String.format("%.1f",inputvals.ymarginratio)+"\n");	
    		output.write("finalplottype : "+ cmbxFinalPlotType.getSelectedIndex()+" : "+inputvals.finalplottype+"\n");			
            output.close();
        } catch ( IOException ioerr ) {
             ioerr.printStackTrace();
        }        
	}

	protected void getStoredSettings(){
		// because we are setting screen values directly, this is a tedious mess to set up JRM	
	    try 
	    {
	    	Scanner scanner =  new Scanner(new File("drawtreeInit.txt"));
	        while (scanner.hasNextLine()){
	        	Scanner linescan =  new Scanner( scanner.nextLine());
	        	linescan.useDelimiter(" : ");
	        	String label = linescan.next();
	        	String value = linescan.next();
	        	if ("intree".equals(label)){
	        		txtInputTree.setText(value);
	        	}
	     		else if ("plotfile".equals(label)){
	     			txtPlot.setText(value);
	    		}
	     		else if ("usefont".equals(label)){
	     			cmbxPlotFont.setSelectedIndex(Integer.parseInt(value));
	     		}
	    		else if ("usebranchlengths".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnUseLenY.setSelected(true);
	    				rbtnUseLenN.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnUseLenY.setSelected(false);
	    				rbtnUseLenN.setSelected(true);
	    			}
	    		}	
	    		else if ("labeldirec".equals(label)){
	    			cmbxLabelAngle.setSelectedIndex(Integer.parseInt(value));
	    		}
	            else if ("labelangle".equals(label)){
	            	cmbxAngle.setSelectedIndex(Integer.parseInt(value));
	            }
	            else if ("treegrows".equals(label)){
	            	if ("horizontal".equals(value))
	            	{
	            		rbtnTreeH.setSelected(true);
	            		rbtnTreeV.setSelected(false);
	            	}
	            	else
	            	{
	            		rbtnTreeH.setSelected(false);
	            		rbtnTreeV.setSelected(true);
	            	}            		
	            }
	    		else if ("treerotation".equals(label)){
	    			txtTreeRotation.setText(value);
	    		}
	    		else if ("treearc".equals(label)){
	    			txtTreeArc.setText(value);
	    		}
	    		else if ("iterationkind".equals(label)){
	    			cmbxIterate.setSelectedIndex(Integer.parseInt(value));
	    		}	
	    		else if ("iterationcount".equals(label)){
	    			txtIterationCount.setText(value);
	    		}
	    		else if ("regularizeangles".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnRegAngleY.setSelected(true);
	    				rbtnRegAngleN.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnRegAngleY.setSelected(false);
	    				rbtnRegAngleN.setSelected(true);
	    			}
	    		}
	    		else if ("avoidlabeloverlap".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnAvoidOverY.setSelected(true);
	    				rbtnAvoidOverN.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnAvoidOverY.setSelected(false);
	    				rbtnAvoidOverN.setSelected(true);
	    			}
	    		}
	    		else if ("branchrescale".equals(label)){
	    			cmbxRescale.setSelectedIndex(Integer.parseInt(value));
	    		}
	    		else if ("branchscaler".equals(label)){
	    			txtBranchScale.setText(value);
	    		}
	    		else if ("relcharhgt".equals(label)){
	    			txtRelCharHgt.setText(value);
	    		}
	    		else if ("xmarginratio".equals(label)){
	    			txtXMarginRatio.setText(value);
	    		}
	    		else if ("ymarginratio".equals(label)){
	    			txtYMarginRatio.setText(value);
	    		}
	    		else if ("finalplottype".equals(label)){
	    			cmbxFinalPlotType.setSelectedIndex(Integer.parseInt(value));
	    		}			
	    		else {
	    			String msg = "Unknown label: ";
	    			msg += label;
	    			msg += " with value: ";
	    			msg += value;
	    			msg += " found in drawtreeInit.txt.";
	    			JOptionPane.showMessageDialog(null, msg, "Warning", JOptionPane.WARNING_MESSAGE);			
	    		}
		    }
	    }
		catch (FileNotFoundException e)
		{
			// if it's not there, use the defaults
			resetDefaults();
		} 	
	}

	protected void resetDefaults()
	{
		// reset DrawTree to default values
		txtInputTree.setText("intree");
		txtPlot.setText("plotfile.ps");
		cmbxPlotFont.setSelectedIndex(0);
		rbtnTreeH.setSelected(false);
		rbtnTreeV.setSelected(true);
		rbtnUseLenY.setSelected(true);
		rbtnUseLenN.setSelected(false);
		cmbxLabelAngle.setSelectedIndex(0);
		lblFixedAngleOf.setEnabled(false);
		cmbxAngle.setEnabled(false);
		cmbxAngle.setSelectedIndex(0);
		txtTreeRotation.setText("90.0");
		txtTreeArc.setText("360.0");
		cmbxIterate.setSelectedIndex(0);
		txtIterationCount.setText("100");
		rbtnRegAngleY.setSelected(false);
		rbtnRegAngleN.setSelected(true);
		lblAvoidOverlap.setEnabled(false);
		rbtnAvoidOverY.setSelected(false);
		rbtnAvoidOverY.setEnabled(false);
		rbtnAvoidOverN.setSelected(true);
		rbtnAvoidOverN.setEnabled(false);
		cmbxRescale.setSelectedIndex(0);
		txtRelCharHgt.setText("0.3333");
		lblBranchScale.setEnabled(false);
		txtBranchScale.setEnabled(false);
		txtBranchScale.setText("1.0");
		cmbxFinalPlotType.setSelectedIndex(0);
		txtXMarginRatio.setText("0.08");
		txtYMarginRatio.setText("0.08");
	}
}
