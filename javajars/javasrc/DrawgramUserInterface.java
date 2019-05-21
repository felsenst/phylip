package phylip;

import java.awt.EventQueue;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JRadioButton;
import javax.swing.JButton;
import javax.swing.JTextField;
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

public class DrawgramUserInterface {

	public class DrawgramData{
		String  intree;
		String  usefont;
		String  plotfile;
		String  plotfileopt;
		String  treegrows;
		String  treestyle;
		boolean usebranchlengths;
		Double  labelangle;
		boolean scalebranchlength;
		Double  branchlength;
		Double  breadthdepthratio;
		Double  stemltreedratio;
		Double  chhttipspratio;
		String  ancnodes;
		Double  xmarginratio;
		Double  ymarginratio;
		String  librarypath;
		boolean doplot; // false = do preview
		String  finalplottype;
	}

	public enum LastPage{COUNT, SIZE, OVERLAP}
	private String filedir;
	private boolean phylipCall;
	private boolean plotCall;
	private String calledBy;
	
	private String ancNodesCBdefault = new String("Weighted");

	private JFrame frmDrawgramControls;
	private JTextField txtLabelAngle;
	private JLabel lblAngleLabels;
	private JTextField txtBranchLen;
	private JTextField txtDepthBreadth;
	private JTextField txtStemLenTreeDpth;
	private JTextField txtCharHgtTipSp;
	private JRadioButton btnTreeH;
	private JRadioButton btnTreeV;
	private JRadioButton btnUseLenY;
	private JRadioButton btnUseLenN;
	private JRadioButton btnBranchScaleAuto;
	private JLabel txtBranchScale;
	private JLabel lblCm;
	private JComboBox cmbxTreeStyle;
	private JComboBox cmbxAncNodes;
	private JButton btnInputTree;
	private JTextField txtInputTree;
	private JTextField txtPlot;
	private JButton btnPlot;
	private JComboBox cmbxPlotFont;
	private JLabel lblFinalPlotType;
	private JComboBox cmbxFinalPlotType;
	private JButton btnPreview;
	private JButton btnQuit;
	private JButton btnPlotFile;
	private JLabel lblPSFont;
	private JLabel lblTreeGrows;
	private JLabel lblTreeStyle;
	private JLabel lblUseBranchLengths;
	private JLabel lblBranchLenScaling;
	private JLabel lblDepthBreadthOfTree;
	private JLabel lblStemLengthTreeDepth;
	private JLabel lblAncestralNodes;
	private JLabel lblIfPresent;
	private JLabel lblCharHgtTipSpace;
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
	 * Launch tnRead;
	
	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DrawgramUserInterface window = new DrawgramUserInterface(args);
					window.frmDrawgramControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public DrawgramUserInterface(String[] args) {
		phylipCall = false;
		plotCall = false;
		if (args.length > 0)
		{
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

	// event handlers
	protected void TreeGrowToggle(boolean ishoriz) {		
		if (ishoriz){
			btnTreeH.setSelected(true);
			btnTreeV.setSelected(false);
		}	
		else{
			btnTreeH.setSelected(false);
			btnTreeV.setSelected(true);
		}		
	}
	
	protected void BranchLengthToggle(boolean uselength) {		
		if (uselength){
			btnUseLenY.setSelected(true);
			btnUseLenN.setSelected(false);
			cmbxAncNodes.setEnabled(true);
			for (int i=0; i<cmbxAncNodes.getItemCount(); i++)
			{
				if(cmbxAncNodes.getItemAt(i).toString().contains(ancNodesCBdefault))
				{
					cmbxAncNodes.setSelectedIndex(i);
				}
			}
		}	
		else{
			btnUseLenY.setSelected(false);
			btnUseLenN.setSelected(true);
			cmbxAncNodes.setEnabled(false);
			ancNodesCBdefault = cmbxAncNodes.getSelectedItem().toString();
			for (int i=0; i<cmbxAncNodes.getItemCount(); i++)
			{
				if(cmbxAncNodes.getItemAt(i).toString().contains("Centered"))
				{
					cmbxAncNodes.setSelectedIndex(i);
				}
			}
		}		
	}
	
	protected void ScaleAutoToggle(boolean isauto) {		
		if (isauto){
			btnBranchScaleAuto.setSelected(true);
			txtBranchScale.setEnabled(false);
			txtBranchLen.setEnabled(false);
			lblCm.setEnabled(false);
		}	
		else{
			btnBranchScaleAuto.setSelected(false);
			txtBranchScale.setEnabled(true);
			txtBranchLen.setEnabled(true);
			lblCm.setEnabled(true);
		}		
	}
	
	protected void ChooseFile(JTextField file) {		
		JFileChooser fileChooser = new JFileChooser( filedir);

		int option = fileChooser.showOpenDialog(frmDrawgramControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}	
	}
	
	protected void LabelAngleToggle() {		
		if ((cmbxTreeStyle.getSelectedItem().toString()).contains("Circular")){
			lblAngleLabels.setEnabled(false);
			txtLabelAngle.setEnabled(false);
		}	
		else{
			lblAngleLabels.setEnabled(true);
			txtLabelAngle.setEnabled(true);
		}		
	}
	
	protected boolean LaunchDrawgramInterface(Boolean doPlot){
		DrawgramData inputdata = getInputVals();
		inputdata.doplot = doPlot;
		DrawgramInterface dg = new DrawgramInterface();
		return (dg.DrawgramRun(inputdata));
		
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
		
		frmDrawgramControls = new JFrame();
		frmDrawgramControls.getContentPane().setBackground(new Color(204, 255, 255));
		frmDrawgramControls.setBackground(new Color(204, 255, 255));
		frmDrawgramControls.setTitle("Drawgram");
		frmDrawgramControls.setFont(new Font("Arial", Font.BOLD, 13));
		frmDrawgramControls.setBounds(100, 100, 620, 530);
		frmDrawgramControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDrawgramControls.setPreferredSize(new Dimension(frmDrawgramControls.getBounds().width, frmDrawgramControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmDrawgramControls.getPreferredSize());
		frmDrawgramControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmDrawgramControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][10.00,grow][pref!,grow]", "[][][][][][][][][][][][][][]"));
		
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
		txtInputTree.setBounds(131, 10, 333, 28);
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
		txtPlot.setBounds(131, 39, 333, 28);
		panel.add(txtPlot, "cell 1 1 2 1,growx");
		
		lblPSFont = new JLabel("Postscript Font:");
		lblPSFont.setFont(new Font("Arial", Font.BOLD, 13));
		lblPSFont.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPSFont, "flowx,cell 1 2,alignx right");
		
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
	
		lblTreeGrows = new JLabel("Tree grows:");
		lblTreeGrows.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTreeGrows, "flowx,cell 1 3 1 1,alignx right");
		
		btnTreeH = new JRadioButton("Horizontally");
		btnTreeH.setFont(new Font("Arial", Font.BOLD, 13));
		btnTreeH.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				TreeGrowToggle(true);
			}
		});
		panel.add(btnTreeH, "cell 2 3");
		
		btnTreeV = new JRadioButton("Vertically");
		btnTreeV.setFont(new Font("Arial", Font.BOLD, 13));
		btnTreeV.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				TreeGrowToggle(false);
			}
		});
		panel.add(btnTreeV, "cell 2 3");
		
		lblTreeStyle = new JLabel("Tree style:");
		lblTreeStyle.setFont(new Font("Arial", Font.BOLD, 13));
		lblTreeStyle.setHorizontalAlignment(SwingConstants.TRAILING);
		panel.add(lblTreeStyle, "flowx,cell 1 4 1 1,alignx right");
		
		cmbxTreeStyle = new JComboBox();
		cmbxTreeStyle.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxTreeStyle.setModel(new DefaultComboBoxModel(new String[] {"Phenogram", "Cladogram", "Curvogram", "Eurogram", "Swoopogram", "Circular tree"}));
		cmbxTreeStyle.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				LabelAngleToggle();
			}
		});
		panel.add(cmbxTreeStyle, "cell 2 4,growx");
		
		lblUseBranchLengths = new JLabel("Use branch lengths:");
		lblUseBranchLengths.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseBranchLengths, "flowx,cell 1 5 1 1,alignx right");
		
		btnUseLenY = new JRadioButton("Yes");
		btnUseLenY.setFont(new Font("Arial", Font.BOLD, 13));
		btnUseLenY.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				BranchLengthToggle(true);
			}
		});
		panel.add(btnUseLenY, "cell 2 5");
		
		btnUseLenN = new JRadioButton("No");
		btnUseLenN.setFont(new Font("Arial", Font.BOLD, 13));
		btnUseLenN.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				BranchLengthToggle(false);
			}
		});
		panel.add(btnUseLenN, "cell 2 5");
		
		lblIfPresent = new JLabel(" (if present)");
		lblIfPresent.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblIfPresent, "cell 2 5");
		
		lblAngleLabels = new JLabel("Angle of labels:");
		lblAngleLabels.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAngleLabels, "flowx,cell 1 6 1 1,alignx right");
		
		txtLabelAngle = new JTextField();
		txtLabelAngle.setFont(new Font("Arial", Font.PLAIN, 13));
		txtLabelAngle.setColumns(6);
		panel.add(txtLabelAngle, "cell 2 6");
		
		lblBranchLenScaling = new JLabel("Branch length scaling:");
		lblBranchLenScaling.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBranchLenScaling, "flowx,cell 1 7 1 1,alignx right");
		
		btnBranchScaleAuto = new JRadioButton("Automatic");
		btnBranchScaleAuto.setFont(new Font("Arial", Font.BOLD, 13));
		btnBranchScaleAuto.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ScaleAutoToggle(btnBranchScaleAuto.isSelected());
			}
		});
		panel.add(btnBranchScaleAuto, "cell 2 7");
	
		txtBranchScale = new JLabel("Scale");
		txtBranchScale.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(txtBranchScale, "cell 2 7");
		
		txtBranchLen = new JTextField();
		txtBranchLen.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBranchLen.setColumns(6);
		panel.add(txtBranchLen, "cell 2 7");
		
		lblCm = new JLabel("cm  ");
		lblCm.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCm, "cell 2 7");
		
		lblDepthBreadthOfTree = new JLabel("Depth/breadth of tree:");
		lblDepthBreadthOfTree.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDepthBreadthOfTree, "flowx,cell 1 8 1 1,alignx right");
		
		txtDepthBreadth = new JTextField();
		txtDepthBreadth.setFont(new Font("Arial", Font.PLAIN, 13));
		txtDepthBreadth.setColumns(6);
		panel.add(txtDepthBreadth, "cell 2 8");
		
		lblStemLengthTreeDepth = new JLabel("Stem length/tree depth:");
		lblStemLengthTreeDepth.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblStemLengthTreeDepth, "flowx,cell 1 9 1 1,alignx right");
		
		txtStemLenTreeDpth = new JTextField();
		txtStemLenTreeDpth.setFont(new Font("Arial", Font.PLAIN, 13));
		txtStemLenTreeDpth.setColumns(6);
		panel.add(txtStemLenTreeDpth, "cell 2 9");
		
		lblCharHgtTipSpace = new JLabel("Character height/tip space:");
		lblCharHgtTipSpace.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCharHgtTipSpace, "flowx,cell 1 10 1 1,alignx right");
 
		txtCharHgtTipSp = new JTextField();
		txtCharHgtTipSp.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCharHgtTipSp.setColumns(6);
		panel.add(txtCharHgtTipSp, "cell 2 10");
		
		lblAncestralNodes = new JLabel("Ancestral nodes:");
		lblAncestralNodes.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAncestralNodes, "flowx,cell 1 11 1 1,alignx right");
		
		cmbxAncNodes = new JComboBox();
		cmbxAncNodes.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxAncNodes.setModel(new DefaultComboBoxModel(new String[] {"Weighted", "Intermediate", "Centered", "Innermost", "V-shaped"}));
		panel.add(cmbxAncNodes, "cell 2 11,growx");
		
		lblFinalPlotType = new JLabel("Final plot file type:");
		lblFinalPlotType.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblFinalPlotType, "flowx,cell 1 12,alignx right");
		
		cmbxFinalPlotType = new JComboBox();
		cmbxFinalPlotType.setModel(new DefaultComboBoxModel(new String[] {"Postscript", "SVG"}));
		cmbxFinalPlotType.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxFinalPlotType, "cell 2 12,growx");
		
		lblMarginRatios = new JLabel("Page margin ratios:");
		lblMarginRatios.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMarginRatios, "flowx,cell 1 13,alignx right");
		
		lblXMarginRatio = new JLabel("X:");
		lblXMarginRatio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblXMarginRatio, "cell 2 13");
		
		txtXMarginRatio = new JTextField();
		txtXMarginRatio.setColumns(5);
		panel.add(txtXMarginRatio, "cell 2 13");
		
		lblYMarginRatio = new JLabel("Y:");
		lblYMarginRatio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblYMarginRatio, "cell 2 13");
		
		txtYMarginRatio = new JTextField();
		txtYMarginRatio.setColumns(5);
		panel.add(txtYMarginRatio, "cell 2 13");
			
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
				boolean retval = LaunchDrawgramInterface(false);
				if (retval)
				{
					String title = "Preview: " + (String)txtPlot.getText();
					String curDir = System.getProperty("user.dir");
					curDir += "/JavaPreview.ps";
					new DrawPreview(title, curDir);
				}
			}
		});	
		panel.add(btnPreview, "cell 2 14");
		
		btnPlotFile = new JButton("Plot");
		btnPlotFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnPlotFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				@SuppressWarnings("unused")
				boolean retval = LaunchDrawgramInterface(true);
			}
		});
		panel.add(btnPlotFile, "cell 2 14");
		
		btnQuit = new JButton("Quit");
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveSettings();
				if((phylipCall) || (plotCall))
				{
					frmDrawgramControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		panel.add(btnQuit, "cell 2 14");	
	}
	
	protected DrawgramData getInputVals()
	{
		DrawgramData inputdata = new DrawgramData();
				
		inputdata.intree = (String)txtInputTree.getText();
		inputdata.usefont =  cmbxPlotFont.getSelectedItem().toString();
		inputdata.plotfile = (String)txtPlot.getText();
		inputdata.plotfileopt = "wb";
		if (btnTreeH.isSelected())
		{
			inputdata.treegrows = "horizontal";
		}
		else
		{
			inputdata.treegrows = "vertical";					
		}
		
		switch (cmbxTreeStyle.getSelectedIndex()){
			case 0: //Phenogram
				inputdata.treestyle = "phenogram";
				break;
			case 1: //Cladogram
				inputdata.treestyle = "cladogram";
				break;
			case 2: //Curvogram
				inputdata.treestyle = "curvogram";
				break;
			case 3: //Eurogram
				inputdata.treestyle = "eurogram";
				break;
			case 4: //Swoopogram"
				inputdata.treestyle = "swoopogram";
				break;
			case 5: //Circular tree
				inputdata.treestyle = "circular";
				break;
			default:
				inputdata.treestyle = "phenogram";
				break;													
		}
		
		inputdata.usebranchlengths = btnUseLenY.isSelected();
		inputdata.labelangle = new Double(txtLabelAngle.getText());
		inputdata.scalebranchlength = btnBranchScaleAuto.isSelected();
		inputdata.branchlength = new Double(txtBranchLen.getText());
		inputdata.breadthdepthratio = new Double(txtDepthBreadth.getText());
		inputdata.stemltreedratio = new Double(txtStemLenTreeDpth.getText());
		inputdata.chhttipspratio = new Double(txtCharHgtTipSp.getText());
		inputdata.xmarginratio = new Double(txtXMarginRatio.getText());
		inputdata.ymarginratio = new Double(txtYMarginRatio.getText());
		
		switch (cmbxAncNodes.getSelectedIndex()){
			case 0: //Weighted
				inputdata.ancnodes = "weighted";
				break;
			case 1: //Intermediate
				inputdata.ancnodes = "intermediate";
				break;
			case 2: //Centered
				inputdata.ancnodes = "centered";
				break;
			case 3: //Innermost
				inputdata.ancnodes = "inner";
				break;
			case 4: //V-shaped
				inputdata.ancnodes = "vshaped";
				break;
			default:
				inputdata.ancnodes = "weighted";
				break;													
		}
		
		switch (cmbxFinalPlotType.getSelectedIndex()){
			case 0: // postscript
				inputdata.finalplottype = "lw";
				break;
			case 1: // Windows Bitmap
				inputdata.finalplottype = "bmp";
				break;
				/*
			case 1: // svg
				inputdata.finalplottype = "svg";
				break;
				*/
			default:
				inputdata.finalplottype = "lw";
				break;	
		}
		
		inputdata.librarypath = System.getProperty("user.dir"); // hardwired - can be made user enterable if need be
		
		return inputdata;
		
	}
	
	protected void saveSettings(){
		DrawgramData inputvals = getInputVals();
		// there must be a better way to format this output, but this works for the prototype JRM
        try {
            BufferedWriter output = new BufferedWriter(new FileWriter("drawgramInit.txt"));
            output.write("intree : "+inputvals.intree+"\n");
            output.write("usefont : "+cmbxPlotFont.getSelectedIndex()+" : "+inputvals.usefont+"\n");
            output.write("plotfile : "+inputvals.plotfile+"\n");
            //output.write("plotfileopt : "+inputvals.plotfileopt+"\n"); // changes between runs
            output.write("treegrows : "+inputvals.treegrows+"\n");
            output.write("treestyle : "+cmbxTreeStyle.getSelectedIndex()+" : "+inputvals.treestyle+"\n");
            output.write("usebranchlengths : "+String.format("%b",inputvals.usebranchlengths)+"\n");
            output.write("labelangle : "+String.format("%.1f",inputvals.labelangle)+"\n");
    		output.write("scalebranchlength : "+String.format("%b",inputvals.scalebranchlength)+"\n");
    		output.write("branchlength : "+String.format("%.1f",inputvals.branchlength)+"\n");
    		output.write("breadthdepthratio : "+String.format("%.2f",inputvals.breadthdepthratio)+"\n");
    		output.write("stemltreedratio : "+String.format("%.2f",inputvals.stemltreedratio)+"\n");
    		output.write("chhttipspratio : "+String.format("%.3f",inputvals.chhttipspratio)+"\n");
    		output.write("ancnodes : "+cmbxAncNodes.getSelectedIndex()+" : "+inputvals.ancnodes+"\n");
    		output.write("xmarginratio : "+String.format("%.2f",inputvals.xmarginratio)+"\n");
    		output.write("ymarginratio : "+String.format("%.2f",inputvals.ymarginratio)+"\n");
    		//output.write("librarypath : "+inputvals.librarypath+"\n"); // currently hardwired
    		//output.write("doplot : "+String.format("%b",inputvals.doplot)+"\n");  // changes between runs
    		output.write("finalplottype : "+cmbxFinalPlotType.getSelectedIndex()+" : "+inputvals.finalplottype+"\n");

            output.close();
        } catch ( IOException ioerr ) {
             ioerr.printStackTrace();
        }        
	}
	
	protected void getStoredSettings(){
		// because we are setting screen values directly, this is a tedious mess to set up JRM
		
	    try 
	    {
	    	Scanner scanner =  new Scanner(new File("drawgramInit.txt"));
	        while (scanner.hasNextLine()){
	        	Scanner linescan =  new Scanner( scanner.nextLine());
	        	linescan.useDelimiter(" : ");
	        	String label = linescan.next();
	        	String value = linescan.next();
	        	if ("intree".equals(label)){
	        		txtInputTree.setText(value);
	        	}
	     		else if ("usefont".equals(label)){
	     			cmbxPlotFont.setSelectedIndex(Integer.parseInt(value));
	     		}
	     		else if ("plotfile".equals(label)){
	     			txtPlot.setText(value);
	    		}
	     		//else if ("plotfileopt".equals(label)){} // not stored
	     		else if ("treegrows".equals(label)){
	    			if ("vertical".equals(value))
	    			{
	    				btnTreeH.setSelected(false);
	    				btnTreeV.setSelected(true);
	    			}
	    			else
	    			{
	    				btnTreeH.setSelected(true);
	    				btnTreeV.setSelected(false);
	    			}
	     		}
	     		else if ("treestyle".equals(label)){
	     			cmbxTreeStyle.setSelectedIndex(Integer.parseInt(value));
	     		}
	     		else if ("usebranchlengths".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				btnUseLenY.setSelected(true);
	    				btnUseLenN.setSelected(false);
	    			}
	    			else
	    			{
	    				btnUseLenY.setSelected(false);
	    				btnUseLenN.setSelected(true);
	    			}
	     		}
	     		else if ("labelangle".equals(label)){
	     			txtLabelAngle.setText(value);
	     		}
	     		else if ("scalebranchlength".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				btnBranchScaleAuto.setSelected(true);
	    				txtBranchScale.setEnabled(false);
	    				txtBranchLen.setEnabled(false);
	    				lblCm.setEnabled(false);
	    			}
	    			else
	    			{
	    				btnBranchScaleAuto.setSelected(false);
	    				txtBranchScale.setEnabled(true);
	    				txtBranchLen.setEnabled(true);
	    				lblCm.setEnabled(true);
	    			}
	     		}
	     		else if ("branchlength".equals(label)){
	     			txtBranchLen.setText(value);
	     		}
	     		else if ("breadthdepthratio".equals(label)){
	     			txtDepthBreadth.setText(value);
	     		}
	     		else if ("stemltreedratio".equals(label)){
	     			txtStemLenTreeDpth.setText(value);
	     		}
	     		else if ("chhttipspratio".equals(label)){
	     			txtCharHgtTipSp.setText(value);	     			
	     		}
	     		else if ("ancnodes".equals(label)){
	     			cmbxAncNodes.setSelectedIndex(Integer.parseInt(value));
	     		}
	     		else if ("xmarginratio".equals(label)){
	     			txtXMarginRatio.setText(value);		
	     		}
	     		else if ("ymarginratio".equals(label)){
	     			txtYMarginRatio.setText(value);
	     		}
	     		//else if ("librarypath".equals(label)){}  // not stored
	     		//else if ("doplot".equals(label)){}  // not stored
	     		else if ("finalplottype".equals(label)){
	     			cmbxFinalPlotType.setSelectedIndex(Integer.parseInt(value));
	     		}
	    		else {
	    			String msg = "Unknown label: ";
	    			msg += label;
	    			msg += " with value: ";
	    			msg += value;
	    			msg += " found in drawgramInit.txt.";
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
		// reset DrawGram to default values
		txtInputTree.setText("intree");
		txtPlot.setText("plotfile.ps");
		cmbxPlotFont.setSelectedIndex(0);
		btnTreeH.setSelected(true);
		btnTreeV.setSelected(false);
		cmbxTreeStyle.setSelectedIndex(0);
		btnUseLenY.setSelected(true);
		btnUseLenN.setSelected(false);
		txtLabelAngle.setText("90.0");
		btnBranchScaleAuto.setSelected(true);
		txtBranchScale.setEnabled(false);
		txtBranchLen.setEnabled(false);
		txtBranchLen.setText("1.0");
		lblCm.setEnabled(false);
		txtDepthBreadth.setText("0.53");
		txtStemLenTreeDpth.setText("0.05");
		txtCharHgtTipSp.setText("0.333");
		cmbxAncNodes.setSelectedIndex(0);
		cmbxFinalPlotType.setSelectedIndex(0);
		txtXMarginRatio.setText("0.08");
		txtYMarginRatio.setText("0.08");
	}
}
