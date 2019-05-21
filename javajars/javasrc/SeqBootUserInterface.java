package phylip;
import java.awt.EventQueue;
import javax.swing.JFileChooser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import javax.swing.JFrame;
import javax.swing.JButton;
import javax.swing.JTextField;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JOptionPane;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import javax.swing.JComboBox;
import javax.swing.SwingUtilities;

import java.awt.Font;
import javax.swing.SwingConstants;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import com.sun.jna.Library;
import com.sun.jna.Native;

import dnaml.DnaMLUserInterface;

import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class SeqBootUserInterface {
   public interface SeqBoot extends Library {
        public void seqboot(
        		String infile,
        		String factorfile,
        		String wgtsfile,
        		String catsfile,
        		String mixfile,
        		String ancfile,
        		String outfile,
        		String outfileopt,
        		String outfactor,
        		String outfactoropt,
        		String outwgts,
        		String outwgtsopt,
        		String outcats,
        		String outcatsopt,
        		String outmix,
        		String outmixopt,
        		String outanc,
        		String outancopt,
        		String DataType,
        		String Method,
        		double SampleFrac,
        		int BlockSize,
        		int Replicates,
        		int Rseed,
        		boolean UseWgts,
        		boolean UseCats,
        		boolean UseFactors,
        		boolean UseMix,
        		boolean UseAnc,
        		boolean HasEnz,
        		boolean AllAlleles,
        		boolean WriteData,
        		boolean InputInterleaved,
        		String OutFmt,
        		String SeqType,
        		boolean PrintData,
        		boolean DotDiff,
        		boolean PrintInd);
   }

	public class SeqBootData {
		String infile;
		String factorfile;
		String wgtsfile;
		String catsfile;
		String mixfile;
		String ancfile;
		String outfile;
		String outfileopt;
		String outfactor;
		String outfactoropt;
		String outwgts;
		String outwgtsopt;
		String outcats;
		String outcatsopt;
		String outmix;
		String outmixopt;
		String outanc;
		String outancopt;
		String DataType;
		String Method;
		double SampleFrac;
		int BlockSize;
		int Replicates;
		int Rseed;
		boolean UseWgts;
		boolean UseCats;
		boolean UseFactors;
		boolean UseMix;
		boolean UseAnc;
		boolean HasEnz;
		boolean AllAlleles;
		boolean WriteData;
		boolean InputInterleaved;
		String OutFmt;
		String SeqType;
		boolean PrintData;
		boolean DotDiff;
		boolean PrintInd;
	}

	private SeqBootData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;
	private boolean bootstrapCall;

	private DefaultComboBoxModel regularModel;
	private DefaultComboBoxModel fractionModel;

	private JFrame frmSeqBootControls;
	private JButton btnInputFile;
	private JTextField txtInputFile;
	private JButton btnOutputFile;
	private JTextField txtOutputFile;
	private JButton btnWeightFile;
	private JTextField txtWeightFile;
	private JTextField txtFactorFile;
	private JButton btnFactorFile;
	private JButton btnOutFact;
	private JTextField txtOutFact;
	private JButton btnCatFile;
	private JTextField txtCatFile;
	private JButton btnExecute;
	private JButton btnQuit;
	private JTextField txtOutWgt;
	private JTextField txtOutCat;
	private JTextField txtMixFile;
	private JButton btnAncFile;
	private JTextField txtAncFile;
	private JButton btnOutMix;
	private JTextField txtOutMix;
	private JButton btnOutAnc;
	private JTextField txtOutAnc;
	private JLabel lblDataType;
	private JComboBox cmbxDataType; 
	private JLabel lblMethod;
	private JComboBox cmbxMethod;
	private JLabel lblSampFrac;
	private JRadioButton rbtnRegular;
	private JRadioButton rbtnFraction;
	private JLabel lblPercent;
	private JTextField txtPercent;
	private JLabel lblRandNum;
	private JTextField txtRandNum;
	private JLabel lblRandNumCmt;
	private JLabel lblBlkSiz;
	private JTextField txtBlkSize;
	private JLabel lblBlkSizeCmt;
	private JLabel lblReps;
	private JTextField txtReps;
	private JLabel lblRdWgts;
	private JRadioButton rbtnRdWgtsYes;
	private JRadioButton rbtnRdWgtsNo;
	private JLabel lblRdSiteCats;
	private JRadioButton rbtnRdSiteCatsYes;
	private JRadioButton rbtnRdSiteCatsNo;
	private JLabel lblWriteWhat;
	private JRadioButton rbtnDataSets;
	private JRadioButton rbtnWeights;
	private JLabel lblInSeq;
	private JRadioButton rbtnInterleaved;
	private JRadioButton rbtnSeqential;
	private JLabel lblPrintData;
	private JRadioButton rbtnPrintDataYes;
	private JRadioButton rbtnPrintDataNo;
	private JLabel lblPrintInd;
	private JRadioButton rbtnPrintIndYes;
	private JRadioButton rbtnPrintIndNo;
	private JLabel lblFactors;
	private JRadioButton rbtnFactorsYes;
	private JRadioButton rbtnFactorsNo;
	private JLabel lblMixFile;
	private JRadioButton rbtnMixFileYes;
	private JRadioButton rbtnMixFileNo;
	private JLabel lblAncFile;
	private JRadioButton rbtnAncFileYes;
	private JRadioButton rbtnAncFileNo;
	private JLabel lblEnzPres;
	private JRadioButton rbtnEnzPresYes;
	private JRadioButton rbtnEnzPresNo;
	private JLabel lblAlleles;
	private JRadioButton rbtnAllelYes;
	private JRadioButton rbtnAllelNo;
	private JButton btnOutCat;
	private JButton btnOutWgt;
	private JButton btnMixFile;
	private JLabel lblOutputFormat;
	private JComboBox cmbxOutputFormat;
	private JLabel lblSeqType;
	private JComboBox cmbxSeqType;
	private JLabel lblDotDiff;
	private JRadioButton rbtnDotDiffYes;
	private JRadioButton rbtnDotDiffNo;
	
	private JScrollPane scrollPane;
	private JPanel panel;
	private JButton btnDefaults;
	private JButton btnStored;

	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					SeqBootUserInterface window = new SeqBootUserInterface(args);
					window.frmSeqBootControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	protected void ChooseFile(JTextField file) {
		// Construct a new file choose whose default path is the path to this
		// executable, which
		// is returned by System.getProperty("user.dir")

		JFileChooser fileChooser = new JFileChooser(filedir);

		int option = fileChooser.showOpenDialog(frmSeqBootControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}

	protected void InputSeqToggle(boolean isInput) {
		if (isInput) {
			rbtnInterleaved.setSelected(true);
			rbtnSeqential.setSelected(false);
		} else {
			rbtnInterleaved.setSelected(false);
			rbtnSeqential.setSelected(true);
		}
	}

	protected void PrintDataToggle(boolean isPrintData) {
		if (isPrintData) {
			rbtnPrintDataYes.setSelected(true);
			rbtnPrintDataNo.setSelected(false);
			lblDotDiff.setEnabled(true);
			rbtnDotDiffYes.setEnabled(true);
			rbtnDotDiffNo.setEnabled(true);
		} else {
			rbtnPrintDataYes.setSelected(false);
			rbtnPrintDataNo.setSelected(true);
			lblDotDiff.setEnabled(false);
			rbtnDotDiffYes.setEnabled(false);
			rbtnDotDiffNo.setEnabled(false);
		}
	}

	protected void DotDiffToggle(boolean isDot) {
		if (isDot) {
			rbtnDotDiffYes.setSelected(true);
			rbtnDotDiffNo.setSelected(false);
		} else {
			rbtnDotDiffYes.setSelected(false);
			rbtnDotDiffNo.setSelected(true);
		}
	}

	protected void PrintIndToggle(boolean isPrintInd) {
		if (isPrintInd) {
			rbtnPrintIndYes.setSelected(true);
			rbtnPrintIndNo.setSelected(false);
		} else {
			rbtnPrintIndYes.setSelected(false);
			rbtnPrintIndNo.setSelected(true);
		}
	}

	protected void WriteKindToggle(boolean isDataSet) {
		if (isDataSet) {
			rbtnDataSets.setSelected(true);
			rbtnWeights.setSelected(false);
			btnOutputFile.setEnabled(true);
			txtOutputFile.setEnabled(true);
			btnOutWgt.setEnabled(false);
			txtOutWgt.setEnabled(false);
			
		} else {
			rbtnDataSets.setSelected(false);
			rbtnWeights.setSelected(true);
			btnOutputFile.setEnabled(false);
			txtOutputFile.setEnabled(false);
			btnOutWgt.setEnabled(true);
			txtOutWgt.setEnabled(true);
		}
	}
	
	protected void ReadSiteCatsToggle(boolean isYes){
		if (isYes)
		{
			rbtnRdSiteCatsYes.setSelected(true);
			rbtnRdSiteCatsNo.setSelected(false);			
			btnCatFile.setEnabled(true);
			txtCatFile.setEnabled(true);
			btnOutCat.setEnabled(true);
			txtOutCat.setEnabled(true);
		} else {
			rbtnRdSiteCatsYes.setSelected(false);
			rbtnRdSiteCatsNo.setSelected(true);
			btnCatFile.setEnabled(false);
			txtCatFile.setEnabled(false);
			btnOutCat.setEnabled(false);
			txtOutCat.setEnabled(false);
		}
	}
	
	protected void ReadWeightsToggle(boolean isYes){
		if (isYes)
		{
			rbtnRdWgtsNo.setSelected(false);
			rbtnRdWgtsYes.setSelected(true);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
		} else {
			rbtnRdWgtsNo.setSelected(true);
			rbtnRdWgtsYes.setSelected(false);
			btnWeightFile.setEnabled(false);
			txtWeightFile.setEnabled(false);
		}
	}
	
	protected void ReadFactorsToggle(boolean isYes){
		if (isYes)
		{
			rbtnFactorsNo.setSelected(false);
			rbtnFactorsYes.setSelected(true);
			btnFactorFile.setEnabled(true);
			txtFactorFile.setEnabled(true);
			btnOutFact.setEnabled(true);
			txtOutFact.setEnabled(true);
		} else {
			rbtnFactorsNo.setSelected(true);
			rbtnFactorsYes.setSelected(false);
			btnFactorFile.setEnabled(false);
			txtFactorFile.setEnabled(false);
			btnOutFact.setEnabled(false);
			txtOutFact.setEnabled(false);
		}
	}
	
	protected void ReadMixtureToggle(boolean isYes){
		if (isYes)
		{
			rbtnMixFileNo.setSelected(false);
			rbtnMixFileYes.setSelected(true);
			btnMixFile.setEnabled(true);
			txtMixFile.setEnabled(true);
			btnOutMix.setEnabled(true);
			txtOutMix.setEnabled(true);
		} else {
			rbtnMixFileNo.setSelected(true);
			rbtnMixFileYes.setSelected(false);
			btnMixFile.setEnabled(false);
			txtMixFile.setEnabled(false);
			btnOutMix.setEnabled(false);
			txtOutMix.setEnabled(false);
		}
	}
	
	protected void ReadAncestorsToggle(boolean isYes){
		if (isYes)
		{
			rbtnAncFileYes.setSelected(true);
			rbtnAncFileNo.setSelected(false);
			btnAncFile.setEnabled(true);
			txtAncFile.setEnabled(true);
			btnOutAnc.setEnabled(true);
			txtOutAnc.setEnabled(true);
		} else {
			rbtnAncFileYes.setSelected(false);
			rbtnAncFileNo.setSelected(true);
			btnAncFile.setEnabled(false);
			txtAncFile.setEnabled(false);
			btnOutAnc.setEnabled(false);
			txtOutAnc.setEnabled(false);
		}
	}
	
	protected void EnzymesToggle(boolean isYes){
		if (isYes)
		{
			rbtnEnzPresYes.setSelected(true);
			rbtnEnzPresNo.setSelected(false);
		} else {
			rbtnEnzPresYes.setSelected(false);
			rbtnEnzPresNo.setSelected(true);
		}
	}
	
	protected void AllelesToggle(boolean isYes){
		if (isYes)
		{
			rbtnAllelYes.setSelected(true);
			rbtnAllelNo.setSelected(false);
		} else {
			rbtnAllelYes.setSelected(false);
			rbtnAllelNo.setSelected(true);
		}
	}

	protected void SamplingToggle(boolean isRegular){
		int MethodIndex = cmbxMethod.getSelectedIndex();
		if (isRegular)
		{
			rbtnRegular.setSelected(true);
			rbtnFraction.setSelected(false);
			lblPercent.setEnabled(false);
			txtPercent.setEnabled(false);
			cmbxMethod.setModel(regularModel);

		} else {
			rbtnRegular.setSelected(false);
			rbtnFraction.setSelected(true);
			lblPercent.setEnabled(true);
			txtPercent.setEnabled(true);
			cmbxMethod.setModel(fractionModel);
		}
		cmbxMethod.setSelectedIndex(MethodIndex);
	}

	protected void DataMethodLogic(int datatype, int method, int outtype){
		
		// method dependent settings
		if ((method == 0) || (method == 1)) // Bootstrap or Jackknife
		{
			enableSampling(true);
		}
		else
		{
			enableSampling(false);
		}
		
		if (method == 0) // Bootstrap
		{
			enableBlkSize(true);				
		}
		else
		{
			enableBlkSize(false);
		}
		
		if (method < 5) // everything except rewrite data
		{
			enableReps(true);
			enableRandSeed(true);
		}
		else
		{
			enableReps(false);
			enableRandSeed(false);
		}
		
		if ((method == 0) || (method == 1) || (method == 2)) // Bootstrap or Jackknife or permute species for each character
		{
			enableRdWgts(true);
		}
		else
		{
			enableRdWgts(false);
		}
		
		if ((method == 0) || (method == 1)) // Bootstrap or Jackknife
		{
			enableOutputType(true);
		}
		else
		{
			WriteKindToggle(true);// hardwire to get things set up right
			enableOutputType(false);
		}
		
		if (method < 5) // everything except rewrite data
		{
			enableOutputFormat(false);
			enableSeqType(false);
		}
		else
		{
			enableOutputFormat(true);
			enableSeqType(false);
			if (outtype == 0) // not PHYLIP
			{
				enableSeqType(false);		
			}
			else 
			{
				enableSeqType(true);			
			}
		}
	
		// datatype/method dependent settings
		if (datatype == 0) // Molecular sequences
		{
			if ((method == 0) || (method == 1) || (method == 2)) // Bootstrap or Jackknife or permute species for each character
			{
				enableSiteCats(true);
			}
			else
			{
				enableSiteCats(false);
			}
			
			enableFactorFile(false);
			enableMixFile(false);
			enableAncFile(false);
			enableEnzPres(false);
			enableAlleles(false);	
			enableInputSeq(true);
		} 
		else if (datatype == 1) // Discrete Morphology
		{		
			enableSiteCats(false);
			
			if (method == 4)  // Permute within species
			{
				enableFactorFile(false);
			}
			else
			{
				enableFactorFile(true);
			}
			
			if ((method == 0) || (method == 1) || (method == 2)) // Bootstrap or Jackknife or permute species for each character
			{
				enableMixFile(true);
				enableAncFile(true);
			}
			else
			{
				enableMixFile(false);
				enableAncFile(false);
			}
			
			enableEnzPres(false);
			enableAlleles(false);					
			enableInputSeq(false);
		}
		else if (datatype == 2) // Restriction Sites
		{
			enableSiteCats(false);
			enableFactorFile(false);
			enableMixFile(false);
			enableAncFile(false);
			enableEnzPres(true);
			enableAlleles(false);						
			enableInputSeq(true);
		}
		else // datatype == 3 // Gene Frequencies
		{
			enableSiteCats(false);
			enableFactorFile(false);
			enableMixFile(false);
			enableAncFile(false);
			enableEnzPres(false);
			enableAlleles(true);						
			enableInputSeq(false);
		}
	}
	protected void enableSampling(boolean val)
	{
		lblSampFrac.setEnabled(val);
		rbtnRegular.setEnabled(val);
		rbtnFraction.setEnabled(val);
		lblPercent.setEnabled(val);
		txtPercent.setEnabled(val);
	}
		
	protected void enableBlkSize(boolean val)
	{
		lblBlkSiz.setEnabled(val);
		txtBlkSize.setEnabled(val);
		lblBlkSizeCmt.setEnabled(val);
	}
	
	protected void enableReps(boolean val)
	{
		lblReps.setEnabled(val);
		txtReps.setEnabled(val);
	}
	
	protected void enableRandSeed(boolean val)
	{
		lblRandNum.setEnabled(val);
		txtRandNum.setEnabled(val);
		lblRandNumCmt.setEnabled(val);
	}
	
	protected void enableRdWgts(boolean val)
	{
		lblRdWgts.setEnabled(val);
		rbtnRdWgtsYes.setEnabled(val);
		rbtnRdWgtsNo.setEnabled(val);
		if (val)
		{
			btnWeightFile.setEnabled(rbtnRdWgtsYes.isSelected());
			txtWeightFile.setEnabled(rbtnRdWgtsYes.isSelected());
		}
		else
		{
			btnWeightFile.setEnabled(val);
			txtWeightFile.setEnabled(val);			
		}
	}
	
	protected void enableSiteCats(boolean val)
	{
		lblRdSiteCats.setEnabled(val);
		rbtnRdSiteCatsYes.setEnabled(val);
		rbtnRdSiteCatsNo.setEnabled(val);
		if (val)
		{
			btnCatFile.setEnabled(rbtnRdSiteCatsYes.isSelected());
			txtCatFile.setEnabled(rbtnRdSiteCatsYes.isSelected());
			btnOutCat.setEnabled(rbtnRdSiteCatsYes.isSelected());
			txtOutCat.setEnabled(rbtnRdSiteCatsYes.isSelected());
		}
		else
		{
			btnCatFile.setEnabled(val);
			txtCatFile.setEnabled(val);
			btnOutCat.setEnabled(val);
			txtOutCat.setEnabled(val);
		}
	}
	
	protected void enableFactorFile(boolean val)
	{
		lblFactors.setEnabled(val);
		rbtnFactorsNo.setEnabled(val);
		rbtnFactorsYes.setEnabled(val);
		if (val)
		{
			btnFactorFile.setEnabled(rbtnFactorsYes.isSelected());
			txtFactorFile.setEnabled(rbtnFactorsYes.isSelected());
			btnOutFact.setEnabled(rbtnFactorsYes.isSelected());
			txtOutFact.setEnabled(rbtnFactorsYes.isSelected());
		}
		else
		{
			btnFactorFile.setEnabled(val);
			txtFactorFile.setEnabled(val);
			btnOutFact.setEnabled(val);
			txtOutFact.setEnabled(val);
		}
	}
	
	protected void enableMixFile(boolean val)
	{
		lblMixFile.setEnabled(val);
		rbtnMixFileNo.setEnabled(val);
		rbtnMixFileYes.setEnabled(val);
		if (val)
		{			
			btnMixFile.setEnabled(rbtnMixFileYes.isSelected());
			txtMixFile.setEnabled(rbtnMixFileYes.isSelected());
			btnOutMix.setEnabled(rbtnMixFileYes.isSelected());
			txtOutMix.setEnabled(rbtnMixFileYes.isSelected());
		}
		else
		{
			btnMixFile.setEnabled(val);
			txtMixFile.setEnabled(val);
			btnOutMix.setEnabled(val);
			txtOutMix.setEnabled(val);
		}
	}
	
	protected void enableAncFile(boolean val)
	{
		lblAncFile.setEnabled(val);
		rbtnAncFileYes.setEnabled(val);
		rbtnAncFileNo.setEnabled(val);
		if (val)
		{			
			btnAncFile.setEnabled(rbtnAncFileYes.isSelected());
			txtAncFile.setEnabled(rbtnAncFileYes.isSelected());
			btnOutAnc.setEnabled(rbtnAncFileYes.isSelected());
			txtOutAnc.setEnabled(rbtnAncFileYes.isSelected());
		}
		else
		{
			btnAncFile.setEnabled(val);
			txtAncFile.setEnabled(val);
			btnOutAnc.setEnabled(val);
			txtOutAnc.setEnabled(val);
		}
	}
	
	protected void enableEnzPres(boolean val)
	{
	
		lblEnzPres.setEnabled(val);
		rbtnEnzPresYes.setEnabled(val);
		rbtnEnzPresNo.setEnabled(val);
	}
	
	protected void enableAlleles(boolean val)
	{
	
		lblAlleles.setEnabled(val);
		rbtnAllelYes.setEnabled(val);
		rbtnAllelNo.setEnabled(val);
	}
	
	protected void enableOutputType(boolean val)
	{
		lblWriteWhat.setEnabled(val);
		rbtnDataSets.setEnabled(val);
		rbtnWeights.setEnabled(val);
	}
	
	protected void enableInputSeq(boolean val)
	{
		lblInSeq.setEnabled(val);
		rbtnInterleaved.setEnabled(val);
		rbtnSeqential.setEnabled(val);
	}
	
	protected void enableOutputFormat(boolean val)
	{
		lblOutputFormat.setEnabled(val);
		cmbxOutputFormat.setEnabled(val);		
	}
	
	protected void enableSeqType(boolean val)
	{
		lblSeqType.setEnabled(val);
		cmbxSeqType.setEnabled(val);		
	}
	
	protected void DoInit()
	{
		// reset everything if there is an init file
		getStoredSettings();
	}

	/**
	 * Create the application.
	 */
	public SeqBootUserInterface(String[] args) {
		phylipCall = false;
		bootstrapCall = false;
		if (args.length > 0)
		{
			if (args[0].contains("Phylip"))
			{
				phylipCall = true;
			}
			if (args.length > 1)
			{
				if (args[1].contains("bootstrap"))
				{
					bootstrapCall = true;
				}
			}
		}
		initialize();
		DoInit();
		
		if (bootstrapCall){
  			txtOutputFile.setText("SeqBootOutfile");
  		}
	}
	
	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		filedir = System.getProperty("user.dir");
		regularModel = new DefaultComboBoxModel(new String[] {
				"Bootstrap", 
				"Delete-half jackknife", 
                "Permute species for each character", 
                "Permute character order", 
                "Permute within species", 
                "Rewrite data"});
		
		fractionModel = new DefaultComboBoxModel(new String[] {
				"Partial Bootstrap", 
				"Delete-fraction jackknife", 
                "Permute species for each character", 
                "Permute character order", 
                "Permute within species", 
                "Rewrite data"});
		
		frmSeqBootControls = new JFrame();
		frmSeqBootControls.setBackground(new Color(204, 255, 255));
		frmSeqBootControls.setTitle("Seqboot");
		frmSeqBootControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmSeqBootControls.setBounds(100, 50, 1200, 800);
		frmSeqBootControls.setPreferredSize(new Dimension(frmSeqBootControls.getBounds().width, frmSeqBootControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmSeqBootControls.getPreferredSize());
		frmSeqBootControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmSeqBootControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][0.00,grow][pref!,grow][pref!,grow][pref!,grow]", "[][][][]"));
		
		btnInputFile = new JButton("Input File");
		btnInputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtInputFile);
			}
		});
		panel.add(btnInputFile, "cell 0 0,growx");
		
		txtInputFile = new JTextField();
		txtInputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtInputFile, "cell 1 0 4 1,growx");
		
		btnFactorFile = new JButton("Factor File");
		btnFactorFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtFactorFile);
			}
		});
		btnFactorFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnFactorFile, "cell 0 1,growx");
		
		txtFactorFile = new JTextField();
		txtFactorFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtFactorFile, "cell 1 1 4 1,growx");
	
		btnWeightFile = new JButton("Weights File");
		btnWeightFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightFile);
			}
		});
		btnWeightFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnWeightFile, "cell 0 2,growx");
	
		txtWeightFile = new JTextField();
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtWeightFile, "cell 1 2 4 1,growx");

		btnCatFile = new JButton("Categories File");
		btnCatFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtCatFile);
			}
		});
		btnCatFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnCatFile, "cell 0 3,growx");
		
		txtCatFile = new JTextField();
		txtCatFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtCatFile, "cell 1 3 4 1,growx");
		
		btnMixFile = new JButton("Mixture File");
		btnMixFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnMixFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtMixFile);
			}
		});
		panel.add(btnMixFile, "cell 0 4,growx");
		
		txtMixFile = new JTextField();
		txtMixFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtMixFile, "cell 1 4 4 1,growx");
		
		btnAncFile = new JButton("Ancestor File");
		btnAncFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnAncFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtAncFile);
			}
		});
		panel.add(btnAncFile, "cell 0 5,growx");
		
		txtAncFile = new JTextField();
		txtAncFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtAncFile, "cell 1 5 4 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		panel.add(btnOutputFile, "cell 0 6,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutputFile, "cell 1 6 4 1,growx");
		
		btnOutFact = new JButton("Output Factors");
		btnOutFact.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutFact.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutFact);
			}
		});
		panel.add(btnOutFact, "cell 0 7,growx");
		
		txtOutFact = new JTextField();
		txtOutFact.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutFact, "cell 1 7 4 1,growx");
		
		btnOutWgt = new JButton("Output Weights");
		btnOutWgt.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutWgt.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutWgt);
			}
		});
		panel.add(btnOutWgt, "cell 0 8,growx");
		
		txtOutWgt = new JTextField();
		txtOutWgt.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutWgt, "cell 1 8 4 1,growx");
		
		btnOutCat = new JButton("Output Categories");
		btnOutCat.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutCat.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutCat);
			}
		});
		panel.add(btnOutCat, "cell 0 9,growx");
		
		txtOutCat = new JTextField();
		txtOutCat.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutCat, "cell 1 9 4 1,growx");
		
		btnOutMix = new JButton("Output Mixtures");
		btnOutMix.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutMix.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutMix);
			}
		});
		panel.add(btnOutMix, "cell 0 10,growx");
		
		txtOutMix = new JTextField();
		txtOutMix.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutMix, "cell 1 10 4 1,growx");
		
		btnOutAnc = new JButton("Output Ancestors");
		btnOutAnc.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutAnc.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutAnc);
			}
		});
		panel.add(btnOutAnc, "cell 0 11,growx");
		
		txtOutAnc = new JTextField();
		txtOutAnc.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtOutAnc, "cell 1 11 4 1,growx");
		
		lblDataType = new JLabel("Input Data Type:");
		lblDataType.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDataType.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDataType, "flowx,cell 0 12 2 1,alignx right");
		
		cmbxDataType = new JComboBox();
		cmbxDataType.setModel(new DefaultComboBoxModel(new String[] {"Molecular sequences", "Discrete Morphology", "Restriction Sites", "Gene Frequencies"}));
		cmbxDataType.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DataMethodLogic(cmbxDataType.getSelectedIndex(), cmbxMethod.getSelectedIndex(), cmbxOutputFormat.getSelectedIndex());
			}
		});
		panel.add(cmbxDataType, "cell 2 12,growx");
		
		lblMethod = new JLabel("Method:");
		lblMethod.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMethod.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMethod, "flowx,cell 0 13 2 1,alignx right");
		
		cmbxMethod = new JComboBox();
		cmbxMethod.setModel(regularModel);
		cmbxMethod.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DataMethodLogic(cmbxDataType.getSelectedIndex(), cmbxMethod.getSelectedIndex(), cmbxOutputFormat.getSelectedIndex());
			}
		});
		panel.add(cmbxMethod, "cell 2 13,growx");
		
		lblSampFrac = new JLabel("Sampling:");
		lblSampFrac.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSampFrac.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSampFrac, "flowx,cell 0 14 2 1,alignx right");
	
		rbtnRegular = new JRadioButton("Regular");
		rbtnRegular.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRegular.setBackground(new Color(204, 255, 255));
		rbtnRegular.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SamplingToggle(true);
			}
		});
		panel.add(rbtnRegular, "cell 2 14");
		
		rbtnFraction = new JRadioButton("Fraction");
		rbtnFraction.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnFraction.setBackground(new Color(204, 255, 255));
		rbtnFraction.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SamplingToggle(false);
			}
		});
		panel.add(rbtnFraction, "cell 2 14");
		
		txtPercent = new JTextField();
		txtPercent.setFont(new Font("Arial", Font.PLAIN, 13));
		txtPercent.setColumns(4);
		panel.add(txtPercent, "cell 2 14");
		
		lblPercent = new JLabel("%");
		lblPercent.setFont(new Font("Arial", Font.BOLD, 13));
		lblPercent.setBackground(new Color(153, 255, 255));
		panel.add(lblPercent, "cell 2 14");
				
		lblBlkSiz = new JLabel("Block Size:");
		lblBlkSiz.setHorizontalAlignment(SwingConstants.RIGHT);
		lblBlkSiz.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBlkSiz, "flowx,cell 0 15 2 1,alignx right");
		
		txtBlkSize = new JTextField();
		txtBlkSize.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBlkSize.setColumns(6);
		txtBlkSize.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				int blksize = Integer.parseInt(txtBlkSize.getText());
				if (blksize == 1)
				{
					lblBlkSizeCmt.setVisible(true);
				}
				else
				{
					lblBlkSizeCmt.setVisible(false);
				}
			}
		});
		panel.add(txtBlkSize, "cell 2 15");
		
		lblBlkSizeCmt = new JLabel("(regular bootstrap)");
		lblBlkSizeCmt.setHorizontalAlignment(SwingConstants.LEFT);
		lblBlkSizeCmt.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBlkSizeCmt, "cell 2 15");
		
		lblReps = new JLabel("Replicates:");
		lblReps.setHorizontalAlignment(SwingConstants.RIGHT);
		lblReps.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblReps, "flowx,cell 0 16 2 1,alignx right");
		
		txtReps = new JTextField();
		txtReps.setFont(new Font("Arial", Font.PLAIN, 13));
		txtReps.setColumns(6);
		panel.add(txtReps, "cell 2 16");
		
		lblRandNum = new JLabel("Random number seed:");
		lblRandNum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRandNum.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandNum, "flowx,cell 0 17 2 1,alignx right");
		
		txtRandNum = new JTextField();
		txtRandNum.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRandNum.setColumns(6);
		panel.add(txtRandNum, "cell 2 17");
		
		lblRandNumCmt = new JLabel("(must be odd)");
		lblRandNumCmt.setHorizontalAlignment(SwingConstants.LEFT);
		lblRandNumCmt.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRandNumCmt, "cell 2 17");
		
		lblEnzPres = new JLabel("Number of enzymes in input file:");
		lblEnzPres.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblEnzPres, "flowx,cell 0 18 2 1,alignx right");
		
		rbtnEnzPresYes = new JRadioButton("Yes");
		rbtnEnzPresYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnEnzPresYes.setBackground(new Color(204, 255, 255));
		rbtnEnzPresYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EnzymesToggle(true);
			}
		});
		panel.add(rbtnEnzPresYes, "cell 2 18");
		
		rbtnEnzPresNo = new JRadioButton("No");
		rbtnEnzPresNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnEnzPresNo.setBackground(new Color(204, 255, 255));
		rbtnEnzPresNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EnzymesToggle(false);
			}
		});
		panel.add(rbtnEnzPresNo, "cell 2 18");
		
		lblAlleles = new JLabel("All alleles present at each locus:");
		lblAlleles.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAlleles.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAlleles, "flowx,cell 0 19 2 1,alignx right");
		
		rbtnAllelYes = new JRadioButton("Yes");
		rbtnAllelYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAllelYes.setBackground(new Color(204, 255, 255));
		rbtnAllelYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AllelesToggle(true);
			}
		});
		panel.add(rbtnAllelYes, "cell 2 19");
		
		rbtnAllelNo = new JRadioButton("No, one absent at each locus");
		rbtnAllelNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAllelNo.setBackground(new Color(204, 255, 255));
		rbtnAllelNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				AllelesToggle(false);
			}
		});
		panel.add(rbtnAllelNo, "cell 2 19");
		
		lblRdWgts = new JLabel("Read character weights:");
		lblRdWgts.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRdWgts.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRdWgts, "flowx,cell 0 20 2 1,alignx right");
		
		rbtnRdWgtsYes = new JRadioButton("Yes");
		rbtnRdWgtsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRdWgtsYes.setBackground(new Color(204, 255, 255));
		rbtnRdWgtsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadWeightsToggle(true);
			}
		});
		panel.add(rbtnRdWgtsYes, "cell 2 20");
		
		rbtnRdWgtsNo = new JRadioButton("No");
		rbtnRdWgtsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRdWgtsNo.setBackground(new Color(204, 255, 255));
		rbtnRdWgtsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadWeightsToggle(false);
			}
		});
		panel.add(rbtnRdWgtsNo, "cell 2 20");
		
		lblRdSiteCats = new JLabel("Read site categories:");
		lblRdSiteCats.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRdSiteCats.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRdSiteCats, "flowx,cell 0 21 2 1,alignx right");
		
		rbtnRdSiteCatsYes = new JRadioButton("Yes");
		rbtnRdSiteCatsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRdSiteCatsYes.setBackground(new Color(204, 255, 255));
		rbtnRdSiteCatsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadSiteCatsToggle(true);
			}
		});
		panel.add(rbtnRdSiteCatsYes, "cell 2 21");
		
		rbtnRdSiteCatsNo = new JRadioButton("No");
		rbtnRdSiteCatsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRdSiteCatsNo.setBackground(new Color(204, 255, 255));
		rbtnRdSiteCatsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadSiteCatsToggle(false);
			}
		});
		panel.add(rbtnRdSiteCatsNo, "cell 2 21");

		/**** start column 2 ****/
		
		lblAncFile = new JLabel("Read ancestor file:");
		lblAncFile.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAncFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAncFile, "cell 3 12, alignx right");
		
		rbtnAncFileYes = new JRadioButton("Yes");
		rbtnAncFileYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAncFileYes.setBackground(new Color(204, 255, 255));
		rbtnAncFileYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadAncestorsToggle(true);
			}
		});
		panel.add(rbtnAncFileYes, "cell 4 12");
		
		rbtnAncFileNo = new JRadioButton("No");
		rbtnAncFileNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnAncFileNo.setBackground(new Color(204, 255, 255));
		rbtnAncFileNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadAncestorsToggle(false);
			}
		});
		panel.add(rbtnAncFileNo, "cell 4 12");
		
		lblFactors = new JLabel("Use factors information:");
		lblFactors.setHorizontalAlignment(SwingConstants.RIGHT);
		lblFactors.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblFactors, "cell 3 13, alignx right");
		
		rbtnFactorsYes = new JRadioButton("Yes");
		rbtnFactorsYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnFactorsYes.setBackground(new Color(204, 255, 255));
		rbtnFactorsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadFactorsToggle(true);
			}
		});
		panel.add(rbtnFactorsYes, "cell 4 13");
		
		rbtnFactorsNo = new JRadioButton("No");
		rbtnFactorsNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnFactorsNo.setBackground(new Color(204, 255, 255));
		rbtnFactorsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadFactorsToggle(false);
			}
		});
		panel.add(rbtnFactorsNo, "cell 4 13");
		
		lblMixFile = new JLabel("Read mixture file:");
		lblMixFile.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMixFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMixFile, "cell 3 14, alignx right");
		
		rbtnMixFileYes = new JRadioButton("Yes");
		rbtnMixFileYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMixFileYes.setBackground(new Color(204, 255, 255));
		rbtnMixFileYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadMixtureToggle(true);
			}
		});
		panel.add(rbtnMixFileYes, "cell 4 14");
		
		rbtnMixFileNo = new JRadioButton("No");
		rbtnMixFileNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMixFileNo.setBackground(new Color(204, 255, 255));
		rbtnMixFileNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ReadMixtureToggle(false);
			}
		});
		panel.add(rbtnMixFileNo, "cell 4 14");
		
		lblWriteWhat = new JLabel("Write out:");
		lblWriteWhat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblWriteWhat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblWriteWhat, "cell 3 15, alignx right");
		
		rbtnDataSets = new JRadioButton("Data Sets");
		rbtnDataSets.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnDataSets.setBackground(new Color(204, 255, 255));
		rbtnDataSets.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteKindToggle(true);
			}
		});
		panel.add(rbtnDataSets, "cell 4 15");
		
		rbtnWeights = new JRadioButton("Just Weights");
		rbtnWeights.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnWeights.setBackground(new Color(204, 255, 255));
		rbtnWeights.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				WriteKindToggle(false);
			}
		});
		panel.add(rbtnWeights, "cell 4 15");
		
		lblInSeq = new JLabel("Input sequences:");
		lblInSeq.setFont(new Font("Arial", Font.BOLD, 13));
		lblInSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblInSeq, "cell 3 16, alignx right");

		rbtnInterleaved = new JRadioButton("Interleaved");
		rbtnInterleaved.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnInterleaved.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnInterleaved.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		panel.add(rbtnInterleaved, "cell 4 16");

		rbtnSeqential = new JRadioButton("Sequential");
		rbtnSeqential.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnSeqential.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnSeqential.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		panel.add(rbtnSeqential, "cell 4 16");
			
		lblOutputFormat = new JLabel("Output format:");
		lblOutputFormat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblOutputFormat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblOutputFormat, "cell 3 17, alignx right");
		
		cmbxOutputFormat = new JComboBox();
		cmbxOutputFormat.setModel(new DefaultComboBoxModel(new String[] {"PHYLIP", "NEXUS", "XML"}));
		cmbxOutputFormat.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DataMethodLogic(cmbxDataType.getSelectedIndex(), cmbxMethod.getSelectedIndex(), cmbxOutputFormat.getSelectedIndex());
			}
		});
		panel.add(cmbxOutputFormat, "cell 4 17, growx");
		
		lblSeqType = new JLabel("Sequence type:");
		lblSeqType.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSeqType.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSeqType, "cell 3 18, alignx right");
		
		cmbxSeqType = new JComboBox();
		cmbxSeqType.setModel(new DefaultComboBoxModel(new String[] {"DNA", "RNA", "Protein"}));
		panel.add(cmbxSeqType, "cell 4 18, growx");
	
		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintData, "cell 3 19, alignx right");
	
		rbtnPrintDataYes = new JRadioButton("Yes");
		rbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		panel.add(rbtnPrintDataYes, "cell 4 19");
	
		rbtnPrintDataNo = new JRadioButton("No");
		rbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		panel.add(rbtnPrintDataNo, "cell 4 19");
	
		lblDotDiff = new JLabel("Use dot-differencing to display:");
		lblDotDiff.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDotDiff, "cell 3 20, alignx right");
	
		rbtnDotDiffYes = new JRadioButton("Yes");
		rbtnDotDiffYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnDotDiffYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDotDiffYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(true);
			}
		});
		panel.add(rbtnDotDiffYes, "cell 4 20");
	
		rbtnDotDiffNo = new JRadioButton("No");
		rbtnDotDiffNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnDotDiffNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDotDiffNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(false);
			}
		});
		panel.add(rbtnDotDiffNo, "cell 4 20");
	
		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintInd, "cell 3 21, alignx right");
	
		rbtnPrintIndYes = new JRadioButton("Yes");
		rbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		panel.add(rbtnPrintIndYes, "cell 4 21");
	
		rbtnPrintIndNo = new JRadioButton("No");
		rbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		panel.add(rbtnPrintIndNo, "cell 4 21");
		
		btnDefaults = new JButton("Restore Defaults");
		btnDefaults.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				resetDefaults();
			}
		});
		btnDefaults.setFont(new Font("Arial", Font.BOLD, 13));	
		panel.add(btnDefaults, "cell 0 22");
		
		btnStored = new JButton("Read Init file");
		btnStored.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				getStoredSettings();
			}
		});
		btnStored.setFont(new Font("Arial", Font.BOLD, 13));	
		panel.add(btnStored, "cell 1 22");

		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = getInputVals();

				btnExecute.setEnabled(false);	
				String title = "Seqboot Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread seqBootThread = new Thread() {
						public void run() {
							runSeqBootThreads();
						}
			  	    };
			  	 	seqBootThread.start();
				}
				btnExecute.setEnabled(true);
				
				if(bootstrapCall)
				{
					saveSettings();
					final String[] bootstrap = new String[] {"Phylip", "bootstrap"};
					try {
						DnaMLUserInterface.main(bootstrap);
					} catch (Exception err) {
						err.printStackTrace();
					}						
					frmSeqBootControls.dispose();
				}
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 4 22");
		
		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				saveSettings();
				if(phylipCall)
				{						
					frmSeqBootControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 4 22");
	}
	
	public boolean checkInputVals(){
		
		// check files
		TestFileNames test = new TestFileNames();
		if (inputvals.WriteData)
		{
			if (!test.DuplicateFileNames(inputvals.infile, "Input", inputvals.outfile, "Output"))
			{			
				return false;		
			}
		}
		
		if (inputvals.UseFactors)
		{
			if (!test.DuplicateFileNames(inputvals.factorfile, "Factor", inputvals.outfactor, "Outfactors"))
			{			
				return false;		
			}
		}

		if (inputvals.UseWgts)
		{
			if (!test.DuplicateFileNames(inputvals.wgtsfile, "Weights", inputvals.outwgts, "Outweights"))
			{			
				return false;		
			}
		}
		
		if (inputvals.UseCats)
		{
			if (!test.DuplicateFileNames(inputvals.catsfile, "Catigories", inputvals.outcats, "Outcatigories"))
			{			
				return false;		
			}
		}
		
		if (inputvals.UseMix)
		{
			if (!test.DuplicateFileNames(inputvals.mixfile, "Mixture", inputvals.outmix, "Outmixture"))
			{			
				return false;		
			}
		}
		
		if (inputvals.UseAnc)
		{
			if (!test.DuplicateFileNames(inputvals.ancfile, "Ancestor", inputvals.outanc, "Outancestors"))
			{			
				return false;		
			}
		}
		
		if (!test.FileAvailable(inputvals.infile, "Input"))
		{
			return false;
		}
		
		if (inputvals.UseFactors)
		{
			if (!test.FileAvailable(inputvals.factorfile, "Factors"))
			{
				return false;
			}
		}
		
		if (inputvals.UseWgts)
		{
			if (!test.FileAvailable(inputvals.wgtsfile, "Weights"))
			{
				return false;
			}
		}
		
		if (inputvals.UseCats)
		{
			if (!test.FileAvailable(inputvals.catsfile, "Categories"))
			{
				return false;
			}
		}
		
		if (inputvals.UseMix)
		{
			if (!test.FileAvailable(inputvals.mixfile, "Mixture"))
			{
				return false;
			}
		}
		
		if (inputvals.UseAnc)
		{
			if (!test.FileAvailable(inputvals.ancfile, "Ancestor"))
			{
				return false;
			}
		}		
		
		String opt;
		if(inputvals.WriteData)
		{
			opt = test.FileAlreadyExists(inputvals.outfile, "Outfile");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inputvals.outfileopt = opt;
			}
		}
		else
		{
			opt = test.FileAlreadyExists(inputvals.outwgts, "Outweights");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inputvals.outwgtsopt = opt;
			}
		}
		
		if (inputvals.UseFactors)
		{
			opt = test.FileAlreadyExists(inputvals.outfactor, "Outfactor");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inputvals.outfactoropt = opt;
			}
		}
				
		if (inputvals.UseCats)
		{
			opt = test.FileAlreadyExists(inputvals.outcats, "Outcatigories");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inputvals.outcatsopt = opt;
			}
		}
		
		if (inputvals.UseMix)
		{
			opt = test.FileAlreadyExists(inputvals.outmix, "Outmixture");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inputvals.outmixopt = opt;
			}
		}
		
		if (inputvals.UseAnc)
		{
			opt = test.FileAlreadyExists(inputvals.outanc, "Outancestors");
			if (opt == "q")
			{
				return false;
			}
			else
			{
				inputvals.outancopt = opt;
			}
		}		

		// check data
		if (inputvals.SampleFrac < 0.0) {
			String msg1 = "Input value: Sampling fraction cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.SampleFrac > 100) {
			String msg1 = "Input value:  Sampling fraction  cannot greater than 100.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if ((inputvals.Rseed % 2) == 0)
		{
			String msg1 = "Random number seed must be odd.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;			
		}
		
		if ((inputvals.UseFactors) && (inputvals.Method == "pspecie"))
		{
			String msg1 = "Cannot use factors when permuting within a species.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;			
		}
		
		if (inputvals.DataType != "seqs")
		{
			if(inputvals.OutFmt == "xml")
			{
				String msg1 = "XML output is not available for ";
				msg1 += inputvals.DataType;
				msg1 += " data.";		
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;			
			}
			
			if (inputvals.UseCats)
			{
				String msg1 = "Categories cannot be used with  ";
				msg1 += inputvals.DataType;
				msg1 += " data.";		
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;						
			}
		}
		
		if (inputvals.DataType != "morphology")
		{
			if (inputvals.UseMix)
			{
				String msg1 = "Mixture file cannot be used with  ";
				msg1 += inputvals.DataType;
				msg1 += " data.";		
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;						
			}
			if (inputvals.UseAnc)
			{
				String msg1 = "Ancestors file cannot be used with  ";
				msg1 += inputvals.DataType;
				msg1 += " data.";		
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;						
			}
		}		

		return true;
	}
	protected void runSeqBootThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("seqboot", SeqBoot.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("SeqBoot");
    		return;
    	}
		try 
		{
	  	    Thread seqBootRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			SeqBoot seqboot = (SeqBoot) Native.loadLibrary("seqboot", SeqBoot.class);
		  			seqboot.seqboot(
		  	        		inputvals.infile,
		  	        		inputvals.factorfile,
		  	        		inputvals.wgtsfile,
		  	        		inputvals.catsfile,
		  	        		inputvals.mixfile,
		  	        		inputvals.ancfile,
		  	        		inputvals.outfile,
		  	        		inputvals.outfileopt,
		  	        		inputvals.outfactor,
		  	        		inputvals.outfactoropt,
		  	        		inputvals.outwgts,
		  	        		inputvals.outwgtsopt,
		  	        		inputvals.outcats,
		  	        		inputvals.outcatsopt,
		  	        		inputvals.outmix,
		  	        		inputvals.outmixopt,
		  	        		inputvals.outanc,
		  	        		inputvals.outancopt,
		  	        		inputvals.DataType,
		  	        		inputvals.Method,
		  	        		inputvals.SampleFrac,
		  	        		inputvals.BlockSize,
		  	        		inputvals.Replicates,
		  	        		inputvals.Rseed,
		  	        		inputvals.UseWgts,
		  	        		inputvals.UseCats,
		  	        		inputvals.UseFactors,
		  	        		inputvals.UseMix,
		  	        		inputvals.UseAnc,
		  	        		inputvals.HasEnz,
		  	        		inputvals.AllAlleles,
		  	        		inputvals.WriteData,
		  	        		inputvals.InputInterleaved,
		  	        		inputvals.OutFmt,
		  	        		inputvals.SeqType,
		  	        		inputvals.PrintData,
		  	        		inputvals.DotDiff,
		  	        		inputvals.PrintInd);
			  	    };
	  	    };
	  	    seqBootRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (seqBootRunThread.isAlive());
	  	    }
		} 
		catch (InterruptedException e) 
		{
			if (inputvals.PrintInd)
			{
				updateProgress();
			}
		}
	}
	
	protected SeqBootData getInputVals()
	{
		SeqBootData inputvals = new SeqBootData();
		inputvals.infile = txtInputFile.getText();
		inputvals.factorfile = txtFactorFile.getText();
		inputvals.wgtsfile = txtWeightFile.getText();
		inputvals.catsfile = txtCatFile.getText();
		inputvals.mixfile = txtMixFile.getText();
		inputvals.ancfile = txtAncFile.getText();
		inputvals.outfile = txtOutputFile.getText();
		inputvals.outfileopt = "w";
		inputvals.outfactor = txtOutFact.getText();
		inputvals.outfactoropt = "w";
		inputvals.outwgts = txtOutWgt.getText();
		inputvals.outwgtsopt = "w";
		inputvals.outcats = txtOutCat.getText();
		inputvals.outcatsopt = "w";
		inputvals.outmix = txtOutMix.getText();
		inputvals.outmixopt = "w";
		inputvals.outanc = txtOutAnc.getText();
		inputvals.outancopt = "w";
		if (cmbxDataType.getSelectedIndex() == 0)
		{
			inputvals.DataType = "seqs";
		}
		else if (cmbxDataType.getSelectedIndex() == 1)
		{
			inputvals.DataType = "morphology";
		}
		else if (cmbxDataType.getSelectedIndex() == 2)
		{
			inputvals.DataType = "restsites";
		}
		else //(cmbxDataType.getSelectedIndex() == 3)
		{
			inputvals.DataType = "genefreqs";
		}

		if (cmbxMethod.getSelectedIndex() == 0)
		{
			if(rbtnRegular.isSelected())
			{
				inputvals.Method = "boot";
			}
			else
			{
				inputvals.Method = "partbt";
			}
		}
		else if (cmbxMethod.getSelectedIndex() == 1)
		{
			if(rbtnRegular.isSelected())
			{
				inputvals.Method = "jack";
			}
			else
			{
				inputvals.Method = "partjk";
			}
		
		}
		else if (cmbxMethod.getSelectedIndex() == 2)
		{
			inputvals.Method = "pseachar";
		}
		else if (cmbxMethod.getSelectedIndex() == 3)
		{
			inputvals.Method = "pcharord";
		}
		else if (cmbxMethod.getSelectedIndex() == 4)
		{
			inputvals.Method = "pspecie";
		}
		else //if (cmbxMethod.getSelectedIndex() == 5)
		{
			inputvals.Method = "rewrite";
		}

		inputvals.SampleFrac = Double.parseDouble(txtPercent.getText());
		inputvals.BlockSize = Integer.parseInt(txtBlkSize.getText());
		inputvals.Replicates = Integer.parseInt(txtReps.getText());
		inputvals.Rseed = Integer.parseInt(txtRandNum.getText());
		inputvals.UseWgts = rbtnRdWgtsYes.isSelected();
		inputvals.UseCats = rbtnRdSiteCatsYes.isSelected();
		inputvals.UseFactors = rbtnFactorsYes.isSelected();
		inputvals.UseMix = rbtnMixFileYes.isSelected();
		inputvals.UseAnc = rbtnAncFileYes.isSelected();
		inputvals.HasEnz = rbtnEnzPresYes.isSelected();
		inputvals.AllAlleles = rbtnAllelYes.isSelected();
		inputvals.WriteData = rbtnDataSets.isSelected();
		inputvals.InputInterleaved = rbtnInterleaved.isSelected();
		if (cmbxOutputFormat.getSelectedIndex() == 0)
		{
			inputvals.OutFmt = "phylip";
		}
		else if (cmbxOutputFormat.getSelectedIndex() == 1)
		{
			inputvals.OutFmt = "nexus";
		}
		else //if (cmbxOutputFormat.getSelectedIndex() == 2)
		{
			inputvals.OutFmt = "xml";
		}
		if (cmbxSeqType.getSelectedIndex() == 0)
		{
			inputvals.SeqType = "DNA";
		}
		else if (cmbxSeqType.getSelectedIndex() == 1)
		{
			inputvals.SeqType = "RNA";
		}
		else //if (cmbxSeqType.getSelectedIndex() == 2)
		{
			inputvals.SeqType = "Protein";
		}
		inputvals.PrintData = rbtnPrintDataYes.isSelected();
		inputvals.DotDiff = rbtnDotDiffYes.isSelected();
		inputvals.PrintInd = rbtnPrintIndYes.isSelected();
		
		return inputvals;
	}
	
	protected void saveSettings(){
		inputvals = getInputVals();
		// there must be a better way to format this output, but this works for the prototype JRM
        try {
            BufferedWriter output = new BufferedWriter(new FileWriter("seqbootInit.txt"));
            output.write("infile : "+inputvals.infile+"\n");
            output.write("factorfile : "+inputvals.factorfile+"\n");
            output.write("wgtsfile : "+inputvals.wgtsfile+"\n");
            output.write("catsfile : "+inputvals.catsfile+"\n");
      		output.write("mixfile : "+inputvals.mixfile+"\n");
      		output.write("ancfile : "+inputvals.ancfile+"\n");
      		output.write("outfile : "+inputvals.outfile+"\n");
      		//output.write("outfileopt : "+inputvals.outfileopt+"\n"); //makes no sense to save, can change between runs JRM
      		output.write("outfactor : "+inputvals.outfactor+"\n");
      		//output.write("outfactoropt : "+inputvals.outfactoropt+"\n"); //makes no sense to save, can change between runs JRM
      		output.write("outwgts : "+inputvals.outwgts+"\n");
      		//output.write("outwgtsopt : "+inputvals.outwgtsopt+"\n"); //makes no sense to save, can change between runs JRM
      		output.write("outcats : "+inputvals.outcats+"\n");
      		//output.write("outcatsopt : "+inputvals.outcatsopt+"\n"); //makes no sense to save, can change between runs JRM
      		output.write("outmix : "+inputvals.outmix+"\n");
      		//output.write("outmixopt : "+inputvals.outmixopt+"\n"); //makes no sense to save, can change between runs JRM
      		output.write("outanc : "+inputvals.outanc+"\n");
      		//output.write("outancopt : "+inputvals.outancopt+"\n"); //makes no sense to save, can change between runs JRM
      		output.write("DataType : "+inputvals.DataType+"\n");
      		output.write("Method : "+inputvals.Method+"\n");
      		output.write("SampleFrac : "+inputvals.SampleFrac+"\n");
      		output.write("BlockSize : "+inputvals.BlockSize+"\n");
      		output.write("Replicates : "+inputvals.Replicates+"\n");
      		output.write("Rseed : "+inputvals.Rseed+"\n");
      		output.write("UseWgts : "+String.format("%b",inputvals.UseWgts)+"\n");
      		output.write("UseCats : "+String.format("%b",inputvals.UseCats)+"\n");
      		output.write("UseFactors : "+String.format("%b",inputvals.UseFactors)+"\n");
      		output.write("UseMix : "+String.format("%b",inputvals.UseMix)+"\n");
      		output.write("UseAnc : "+String.format("%b",inputvals.UseAnc)+"\n");
      		output.write("HasEnz : "+String.format("%b",inputvals.HasEnz)+"\n");
      		output.write("AllAlleles : "+String.format("%b",inputvals.AllAlleles)+"\n");
      		output.write("WriteData : "+String.format("%b",inputvals.WriteData)+"\n");
      		output.write("InputInterleaved : "+String.format("%b",inputvals.InputInterleaved)+"\n");
      		output.write("OutFmt : "+inputvals.OutFmt+"\n");
      		output.write("SeqType : "+inputvals.SeqType+"\n");
      		output.write("PrintData : "+String.format("%b",inputvals.PrintData)+"\n");
      		output.write("DotDiff : "+String.format("%b",inputvals.DotDiff)+"\n");
      		output.write("PrintInd : "+String.format("%b",inputvals.PrintInd)+"\n");
            output.close();
        } catch ( IOException ioerr ) {
             ioerr.printStackTrace();
        }        
	}	
	
	protected void getStoredSettings(){
		// because we are setting screen values directly, this is a tedious mess to set up JRM
	    try 
	    {
	    	Scanner scanner =  new Scanner(new File("seqbootInit.txt"));
	        while (scanner.hasNextLine()){
	        	Scanner linescan =  new Scanner( scanner.nextLine());
	        	linescan.useDelimiter(" : ");
	        	String label = linescan.next();
	        	String value = linescan.next();
	        	if ("infile".equals(label)){
	        		txtInputFile.setText(value);
	        	}
	        	else if("factorfile".equals(label)){
	        		txtFactorFile.setText(value);
	        	}
	            else if("wgtsfile".equals(label)){
	            	txtWeightFile.setText(value);
	        	}
	            else if("catsfile".equals(label)){
	            	txtCatFile.setText(value);
	        	}
	      		else if("mixfile".equals(label)){
	      			txtMixFile.setText(value);
	        	}
	      		else if("ancfile".equals(label)){
	      			txtAncFile.setText(value);
	        	}
	      		else if("outfile".equals(label)){
	      			txtOutputFile.setText(value);
	        	}
	      		else if("outfactor".equals(label)){
	      			txtOutFact.setText(value);
	        	}
	      		else if("outwgts".equals(label)){
	      			txtOutWgt.setText(value);
	        	}
	      		else if("outcats".equals(label)){
	      			txtOutCat.setText(value);
	        	}
	      		else if("outmix".equals(label)){
	      			txtOutMix.setText(value);
	        	}
	      		else if("outanc".equals(label)){
	      			txtOutAnc.setText(value);
	        	}
	      		else if("DataType".equals(label)){
					if (value.contains("seqs"))
					{
						cmbxDataType.setSelectedIndex(0);
					}
					else if (value.contains("morphology"))
					{
						cmbxDataType.setSelectedIndex(1);				
					}
					else if (value.contains("restsites"))
					{
						cmbxDataType.setSelectedIndex(2);				
					}
					else //if (value.contains("genefreqs"))
					{
						cmbxDataType.setSelectedIndex(3);				
					}
	        	}
	      		else if("Method".equals(label)){
					if (value.contains("boot"))
					{
						cmbxMethod.setSelectedIndex(0);
						rbtnRegular.setSelected(true);
					}
					else if (value.contains("partbt"))
					{
						cmbxMethod.setSelectedIndex(0);
						rbtnRegular.setSelected(false);
					}
					else if (value.contains("jack"))
					{
						cmbxMethod.setSelectedIndex(1);
						rbtnRegular.setSelected(true);
					}
					else if (value.contains("partjk"))
					{
						cmbxMethod.setSelectedIndex(1);
						rbtnRegular.setSelected(false);
					}
					else if (value.contains("pseachar"))
					{
						cmbxMethod.setSelectedIndex(2);
					}
					else if (value.contains("pcharord"))
					{
						cmbxMethod.setSelectedIndex(3);
					}
					else if (value.contains("pspecie"))
					{
						cmbxMethod.setSelectedIndex(4);
					}
					else //if (value.contains("rewrite"))
					{
						cmbxMethod.setSelectedIndex(5);
					}
	        	}
	      		else if("SampleFrac".equals(label)){
	      			txtPercent.setText(value);
	        	}
	      		else if("BlockSize".equals(label)){
	      			txtBlkSize.setText(value);
	        	}
	      		else if("Replicates".equals(label)){
	      			txtReps.setText(value);
	        	}
	      		else if("Rseed".equals(label)){
	      			txtRandNum.setText(value);
	        	}
	      		else if("UseWgts".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnRdWgtsYes.setSelected(true);
	    				rbtnRdWgtsNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnRdWgtsYes.setSelected(false);
	    				rbtnRdWgtsNo.setSelected(true);
	    			}
	        	}
	      		else if("UseCats".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnRdSiteCatsYes.setSelected(true);
	    				rbtnRdSiteCatsNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnRdSiteCatsYes.setSelected(false);
	    				rbtnRdSiteCatsNo.setSelected(true);
	    			}
	        	}
	      		else if("UseFactors".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnFactorsYes.setSelected(true);
	    				rbtnFactorsNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnFactorsYes.setSelected(false);
	    				rbtnFactorsNo.setSelected(true);
	    			}
	        	}
	      		else if("UseMix".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnMixFileYes.setSelected(true);
	    				rbtnMixFileNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnMixFileYes.setSelected(false);
	    				rbtnMixFileNo.setSelected(true);
	    			}
	        	}
	      		else if("UseAnc".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnAncFileYes.setSelected(true);
	    				rbtnAncFileNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnAncFileYes.setSelected(false);
	    				rbtnAncFileNo.setSelected(true);
	    			}
	        	}
	      		else if("HasEnz".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnEnzPresYes.setSelected(true);
	    				rbtnEnzPresNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnEnzPresYes.setSelected(false);
	    				rbtnEnzPresNo.setSelected(true);
	    			}
	        	}
	      		else if("AllAlleles".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnAllelYes.setSelected(true);
	    				rbtnAllelNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnAllelYes.setSelected(false);
	    				rbtnAllelNo.setSelected(true);
	    			}
	        	}
	      		else if("WriteData".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnDataSets.setSelected(true);
	    			}
	    			else
	    			{
	    				rbtnDataSets.setSelected(false);
	    			}
	        	}
	      		else if("InputInterleaved".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnInterleaved.setSelected(true);
	    			}
	    			else
	    			{
	    				rbtnInterleaved.setSelected(false);
	    			}
	        	}
	      		else if("OutFmt".equals(label)){
					if (value.contains("phylip"))
					{
						cmbxOutputFormat.setSelectedIndex(0);
					}
					else if (value.contains("nexus"))
					{
						cmbxOutputFormat.setSelectedIndex(1);
					}
					else //if (value.contains("xml"))
					{
						cmbxOutputFormat.setSelectedIndex(2);
					}
	        	}
	      		else if("SeqType".equals(label)){
					if (value.contains("DNA"))
					{
						cmbxSeqType.setSelectedIndex(0);
					}
					else if (value.contains("RNA"))
					{
						cmbxSeqType.setSelectedIndex(1);
					}
					else //if (value.contains("Protein"))
					{
						cmbxSeqType.setSelectedIndex(2);
					}
	        	}
	      		else if("PrintData".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnPrintDataYes.setSelected(true);
	    				rbtnPrintDataNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnPrintDataYes.setSelected(false);
	    				rbtnPrintDataNo.setSelected(true);
	    			}
	        	}
	      		else if("DotDiff".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnDotDiffYes.setSelected(true);
	    				rbtnDotDiffNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnDotDiffYes.setSelected(false);
	    				rbtnDotDiffNo.setSelected(true);
	    			}
	        	}
	      		else if("PrintInd".equals(label)){
	    			if ("true".equals(value))
	    			{
	    				rbtnPrintIndYes.setSelected(true);
	    				rbtnPrintIndNo.setSelected(false);
	    			}
	    			else
	    			{
	    				rbtnPrintIndYes.setSelected(false);
	    				rbtnPrintIndNo.setSelected(true);
	    			}
	        	}
	    		else {
	    			String msg = "Unknown label: ";
	    			msg += label;
	    			msg += " with value: ";
	    			msg += value;
	    			msg += " found in seqbootInit.txt.";
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
		// reset SeqBoot defaults
		txtInputFile.setText("infile");
		btnFactorFile.setEnabled(false);
		txtFactorFile.setEnabled(false);
		txtFactorFile.setText("factor");
		btnWeightFile.setEnabled(false);
		txtWeightFile.setEnabled(false);
		txtWeightFile.setText("weightfile");
		btnCatFile.setEnabled(false);
		txtCatFile.setEnabled(false);
		txtCatFile.setText("catfile");
		btnMixFile.setEnabled(false);
		txtMixFile.setEnabled(false);
		txtMixFile.setText("mixturefile");
		btnAncFile.setEnabled(false);
		txtAncFile.setEnabled(false);
		txtAncFile.setText("ancestorfile");
		txtOutputFile.setText("outfile");
		btnOutFact.setEnabled(false);
		txtOutFact.setEnabled(false);
		txtOutFact.setText("outfactors");
		btnOutWgt.setEnabled(false);
		txtOutWgt.setEnabled(false);
		txtOutWgt.setText("outweights");
		btnOutCat.setEnabled(false);
		txtOutCat.setEnabled(false);
		txtOutCat.setText("outcategories");
		btnOutMix.setEnabled(false);
		txtOutMix.setEnabled(false);
		txtOutMix.setText("outmixture");
		btnOutAnc.setEnabled(false);
		txtOutAnc.setEnabled(false);
		txtOutAnc.setText("outancestors");
		cmbxDataType.setSelectedIndex(0);
		cmbxMethod.setSelectedIndex(0);
		rbtnRegular.setSelected(true);
		rbtnFraction.setSelected(false);
		txtPercent.setEnabled(false);
		txtPercent.setText("50");
		lblPercent.setEnabled(false);
		txtBlkSize.setText("1");
		lblBlkSizeCmt.setVisible(true);
		txtReps.setText("100");
		txtRandNum.setText("1");
		lblEnzPres.setEnabled(false);
		rbtnEnzPresYes.setEnabled(false);
		rbtnEnzPresYes.setSelected(false);
		rbtnEnzPresNo.setEnabled(false);
		rbtnEnzPresNo.setSelected(true);
		lblAlleles.setEnabled(false);
		rbtnAllelYes.setEnabled(false);
		rbtnAllelYes.setSelected(false);
		rbtnAllelNo.setEnabled(false);
		rbtnAllelNo.setSelected(true);
		rbtnRdWgtsYes.setSelected(false);
		rbtnRdWgtsNo.setSelected(true);
		rbtnRdSiteCatsYes.setSelected(false);
		rbtnRdSiteCatsNo.setSelected(true);
		lblAncFile.setEnabled(false);
		rbtnAncFileYes.setEnabled(false);
		rbtnAncFileYes.setSelected(false);
		rbtnAncFileNo.setEnabled(false);
		rbtnAncFileNo.setSelected(true);
		lblFactors.setEnabled(false);
		rbtnFactorsYes.setEnabled(false);
		rbtnFactorsYes.setSelected(false);
		rbtnFactorsNo.setEnabled(false);
		rbtnFactorsNo.setSelected(true);
		lblMixFile.setEnabled(false);
		rbtnMixFileYes.setEnabled(false);
		rbtnMixFileYes.setSelected(false);
		rbtnMixFileNo.setEnabled(false);
		rbtnMixFileNo.setSelected(true);
		rbtnDataSets.setSelected(true);
		rbtnWeights.setSelected(false);
		rbtnInterleaved.setSelected(true);
		rbtnSeqential.setSelected(false);
		lblOutputFormat.setEnabled(false);
		cmbxOutputFormat.setEnabled(false);
		cmbxOutputFormat.setSelectedIndex(0);
		lblSeqType.setEnabled(false);
		cmbxSeqType.setSelectedIndex(0);
		cmbxSeqType.setEnabled(false);
		rbtnPrintDataYes.setSelected(false);
		rbtnPrintDataNo.setSelected(true);
		lblDotDiff.setEnabled(false);
		rbtnDotDiffYes.setEnabled(false);
		rbtnDotDiffYes.setSelected(true);
		rbtnDotDiffNo.setEnabled(false);
		rbtnDotDiffNo.setSelected(false);
		rbtnPrintIndYes.setSelected(true);
		rbtnPrintIndNo.setSelected(false);
}
		
	private void updateProgress(){
		SwingUtilities.invokeLater(new Runnable(){
			public void run() {
				DisplayProgress oldp = dp;
				dp = new DisplayProgress(inTitle,inCurdir);
				// doing the dispose in this order keeps the main window from flickering
				if (oldp != null) 
				{
					oldp.dispose();
				}
			}
		});
	}				
}
