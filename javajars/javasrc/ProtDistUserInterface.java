package phylip;

import java.awt.EventQueue;
import javax.swing.JFileChooser;
import java.io.File;
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

import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class ProtDistUserInterface {
	public interface ProtDist extends Library {
		public void protdist(
			String infile,
			String wgtfile,
			String catsfile,
			String outfile,
			String outfileopt,
			String Model,
			String GammaDist,
			double FracInvar,
			double CoeffVar,
			boolean OneSite,
			int NumSites,
			
			// these are explicitly named because JNA doesn't pass arrays gracefully
			double SiteRate1,  
			double SiteRate2,  
			double SiteRate3,  
			double SiteRate4,  
			double SiteRate5,  
			double SiteRate6,  
			double SiteRate7,  
			double SiteRate8,  
			double SiteRate9,  
			//
	
			boolean SitesWeighted,
			String GeneCode,
			String AACat,
			double ProbCat,
			double TTratio,
			boolean useEqualBF,
			double BaseFreqA,
			double BaseFreqC,
			double BaseFreqG,
			double BaseFreqTU,
			boolean MultData,
			boolean MultDSet,
			int NumSets,
			boolean InputSeq,
			boolean PrintData,
			boolean DotDiff,
			boolean PrintInd);
	}

	public class ProtDistData {
		String infile;
		String wgtfile;
		String catsfile;
		String outfile;
		String outfileopt;
		String Model;
		String GammaDist;
		double FracInvar;
		double CoeffVar;
		boolean OneSite;
		int NumSites;
		
		// these are explicitly named because JNA doesn't pass arrays gracefully
		double SiteRate1;  
		double SiteRate2;  
		double SiteRate3;  
		double SiteRate4;  
		double SiteRate5;  
		double SiteRate6;  
		double SiteRate7;  
		double SiteRate8;  
		double SiteRate9;  
		//

		boolean SitesWeighted;
		String GeneCode;
		String AACat;
		double ProbCat;
		double TTratio;
		boolean useEqualBF;
		double BaseFreqA;
		double BaseFreqC;
		double BaseFreqG;
		double BaseFreqTU;
		boolean MultData;
		boolean MultDSet;
		int NumSets;
		boolean InputSeq;
		boolean PrintData;
		boolean DotDiff;
		boolean PrintInd;
	}

	private ProtDistData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;

	private double[] SiteRateValues;
	private int lastsiteratecat;

	private JFrame frmProtDistControls;
	private JButton btnInputFile;
	private JTextField txtInputFile;
	private JButton btnOutputFile;
	private JTextField txtOutputFile;
	private JButton btnWeightFile;
	private JTextField txtWeightFile;
	private JButton btnCatFile;
	private JTextField txtCatFile;
	private JLabel lblModel;
	private JComboBox cmbxModel;
	private JLabel lblGamma;
	private JComboBox cmbxGamma;
	private JLabel lblFracInvar;
	private JTextField txtFracInvar;
	private JLabel lblCoeffVar;
	private JTextField txtCoeffVar;
	private JLabel lblCoeffVarNote;
	private JLabel lblCatSites;
	private JRadioButton rdbtnCatYes;
	private JRadioButton rdbtnCatNo;
	private JLabel lblNumCat;
	private JComboBox cmbxNumCat;
	private JLabel lblSiteRate;
	private JTextField txtSiteRate;
	private JLabel lblCatNum;
	private JComboBox cmbxRateCatNum;
	private JLabel lblPosWeight;
	private JRadioButton rbtnPosWeightYes;
	private JRadioButton rbtnPosWeightNo;
	private JLabel lblGenCode;
	private JComboBox cmbxGenCode;
	private JLabel lblAACat;
	private JComboBox cmbxAACat;
	private JLabel lblAAInfo;
	private JLabel lblAddAA;
	private JLabel lblProbChange;
	private JTextField txtProbChange;
	private JLabel lblProbChangeInfo;
	private JLabel lblTTRatio;
	private JTextField txtTTRatio;
	private JLabel lblEqualBF;
	private JRadioButton rdbtnEqualBFYes;
	private JRadioButton rdbtnEqualBFNo;
	private JLabel lblBaseFreq;
	private JLabel lblBaseFreqA;
	private JTextField txtBaseFreqA;
	private JLabel lblBaseFreqC;
	private JTextField txtBaseFreqC;
	private JLabel lblBaseFreqG;
	private JTextField txtBaseFreqG;
	private JLabel lblBaseFreqTU;
	private JTextField txtBaseFreqTU;
	private JLabel lblAnalyzeMultData;
	private JRadioButton rbtnMultDataYes;
	private JRadioButton rbtnMultDataNo;
	private JLabel lblMultData;
	private JRadioButton rbtnDataSets;
	private JRadioButton rbtnWeights;
	private JLabel lblHowManyData;
	private JTextField txtNumSeqs;
	private JLabel lblInputSeq;
	private JRadioButton rdbtnInputSeqYes;
	private JRadioButton rdbtnInputSeqNo;
	private JLabel lblPrintData;
	private JRadioButton rbtnPrintDataYes;
	private JRadioButton rbtnPrintDataNo;
	private JLabel lblDotDiff;
	private JRadioButton rbtnDotDiffYes;
	private JRadioButton rbtnDotDiffNo;
	private JLabel lblPrintInd;
	private JRadioButton rbtnPrintIndYes;
	private JRadioButton rbtnPrintIndNo;
	private JButton btnExecute;
	private JButton btnQuit;
	
	private JScrollPane scrollPane;
	private JPanel panel;

	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					ProtDistUserInterface window = new ProtDistUserInterface(args);
					window.frmProtDistControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmProtDistControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}
	
	
	protected void ModelToggle(int selected){
		// remember selected is 0 based

		if(selected == 8) // Categories
		{	
			lblGenCode.setEnabled(true);
			cmbxGenCode.setEnabled(true);
			lblAACat.setEnabled(true);
			cmbxAACat.setEnabled(true);
			lblAAInfo.setEnabled(true);
			lblAddAA.setEnabled(true);
			lblProbChange.setEnabled(true);
			txtProbChange.setEnabled(true);
			lblProbChangeInfo.setEnabled(true);
			lblTTRatio.setEnabled(true);
			txtTTRatio.setEnabled(true);
			lblEqualBF.setEnabled(true);
			rdbtnEqualBFYes.setEnabled(true);
			rdbtnEqualBFNo.setEnabled(true);
			lblGamma.setEnabled(true);
			cmbxGamma.setEnabled(true);
			if(cmbxGamma.getSelectedIndex() == 2)
			{
				lblFracInvar.setEnabled(true);
				txtFracInvar.setEnabled(true);
			}
			if (rdbtnEqualBFNo.isSelected())
			{
				lblBaseFreq.setEnabled(true);
				lblBaseFreqA.setEnabled(true);
				txtBaseFreqA.setEnabled(true);
				lblBaseFreqC.setEnabled(true);
				txtBaseFreqC.setEnabled(true);
				lblBaseFreqG.setEnabled(true);
				txtBaseFreqG.setEnabled(true);
				lblBaseFreqTU.setEnabled(true);
				txtBaseFreqTU.setEnabled(true);				
			}
			else
			{
				lblBaseFreq.setEnabled(false);
				lblBaseFreqA.setEnabled(false);
				txtBaseFreqA.setEnabled(false);
				lblBaseFreqC.setEnabled(false);
				txtBaseFreqC.setEnabled(false);
				lblBaseFreqG.setEnabled(false);
				txtBaseFreqG.setEnabled(false);
				lblBaseFreqTU.setEnabled(false);
				txtBaseFreqTU.setEnabled(false);				
			}
		}
		else
		{
			lblGenCode.setEnabled(false);
			cmbxGenCode.setEnabled(false);
			lblAACat.setEnabled(false);
			cmbxAACat.setEnabled(false);
			lblAAInfo.setEnabled(false);
			lblAddAA.setEnabled(false);
			lblProbChange.setEnabled(false);
			txtProbChange.setEnabled(false);
			lblProbChangeInfo.setEnabled(false);
			lblTTRatio.setEnabled(false);
			txtTTRatio.setEnabled(false);
			lblEqualBF.setEnabled(false);
			rdbtnEqualBFYes.setEnabled(false);
			rdbtnEqualBFNo.setEnabled(false);
			lblBaseFreq.setEnabled(false);
			lblBaseFreqA.setEnabled(false);
			txtBaseFreqA.setEnabled(false);
			lblBaseFreqC.setEnabled(false);
			txtBaseFreqC.setEnabled(false);
			lblBaseFreqG.setEnabled(false);
			txtBaseFreqG.setEnabled(false);
			lblBaseFreqTU.setEnabled(false);
			txtBaseFreqTU.setEnabled(false);
			if(selected == 0) //Jones-Taylor-Thornton
			{
				lblGamma.setEnabled(true);
				cmbxGamma.setEnabled(true);
				if(cmbxGamma.getSelectedIndex() == 2)
				{
					lblFracInvar.setEnabled(true);
					txtFracInvar.setEnabled(true);
				}
				if(cmbxGamma.getSelectedIndex() != 0)
				{
					lblCoeffVar.setEnabled(true);
					txtCoeffVar.setEnabled(true);
					lblCoeffVarNote.setEnabled(true);
				}
			}
			else if(selected == 1) // Heinkoff/Tiller PMB
			{
				lblGamma.setEnabled(true);
				cmbxGamma.setEnabled(true);
				if(cmbxGamma.getSelectedIndex() == 2)
				{
					lblFracInvar.setEnabled(true);
					txtFracInvar.setEnabled(true);
				}
				if(cmbxGamma.getSelectedIndex() != 0)
				{
					lblCoeffVar.setEnabled(true);
					txtCoeffVar.setEnabled(true);
					lblCoeffVarNote.setEnabled(true);
				}
			}
			else if(selected == 2) // Dayhoff PAM
			{
				lblGamma.setEnabled(true);
				cmbxGamma.setEnabled(true);
				if(cmbxGamma.getSelectedIndex() == 2)
				{
					lblFracInvar.setEnabled(true);
					txtFracInvar.setEnabled(true);
				}
				if(cmbxGamma.getSelectedIndex() != 0)
				{
					lblCoeffVar.setEnabled(true);
					txtCoeffVar.setEnabled(true);
					lblCoeffVarNote.setEnabled(true);
				}
			}
			else if(selected == 3) // mtRev
			{
				lblGamma.setEnabled(true);
				cmbxGamma.setEnabled(true);
				if(cmbxGamma.getSelectedIndex() == 2)
				{
					lblFracInvar.setEnabled(true);
					txtFracInvar.setEnabled(true);
				}
				if(cmbxGamma.getSelectedIndex() != 0)
				{
					lblCoeffVar.setEnabled(true);
					txtCoeffVar.setEnabled(true);
					lblCoeffVarNote.setEnabled(true);
				}
			}
			else if(selected == 4) // Whelan and Goldman (WAG)
			{
				lblGamma.setEnabled(true);
				cmbxGamma.setEnabled(true);
				if(cmbxGamma.getSelectedIndex() == 2)
				{
					lblFracInvar.setEnabled(true);
					txtFracInvar.setEnabled(true);
				}
				if(cmbxGamma.getSelectedIndex() != 0)
				{
					lblCoeffVar.setEnabled(true);
					txtCoeffVar.setEnabled(true);
					lblCoeffVarNote.setEnabled(true);
				}
			}
			else if(selected == 5) // mtMam
			{
				lblGamma.setEnabled(true);
				cmbxGamma.setEnabled(true);
				if(cmbxGamma.getSelectedIndex() == 2)
				{
					lblFracInvar.setEnabled(true);
					txtFracInvar.setEnabled(true);
				}
				if(cmbxGamma.getSelectedIndex() != 0)
				{
					lblCoeffVar.setEnabled(true);
					txtCoeffVar.setEnabled(true);
					lblCoeffVarNote.setEnabled(true);
				}
			}
			else if(selected == 6) // Kimura protein
			{
				lblGamma.setEnabled(false);
				cmbxGamma.setEnabled(false);
				lblFracInvar.setEnabled(false);
				txtFracInvar.setEnabled(false);
				lblCoeffVar.setEnabled(false);
				txtCoeffVar.setEnabled(false);
				lblCoeffVarNote.setEnabled(false);
			}
			else // if(selected == 7) // Similarity table
			{
				lblGamma.setEnabled(false);
				cmbxGamma.setEnabled(false);
				lblFracInvar.setEnabled(false);
				txtFracInvar.setEnabled(false);
				lblCoeffVar.setEnabled(false);
				txtCoeffVar.setEnabled(false);
				lblCoeffVarNote.setEnabled(false);
			}
		}
	}

	protected void GammaToggle(int selected){
		if(selected == 0) //No
		{
			lblFracInvar.setEnabled(false);
			txtFracInvar.setEnabled(false);
			lblCoeffVar.setEnabled(false);
			txtCoeffVar.setEnabled(false);
			lblCoeffVarNote.setEnabled(false);
		}
		else{
			lblCoeffVar.setEnabled(true);
			txtCoeffVar.setEnabled(true);
			lblCoeffVarNote.setEnabled(true);
			if(selected == 1) // Yes
			{
				lblFracInvar.setEnabled(false);
				txtFracInvar.setEnabled(false);
			}
			else // Gamma + invariant sites
			{
				lblFracInvar.setEnabled(true);
				txtFracInvar.setEnabled(true);
			}
		}
	}
	
	protected void CatToggle(boolean isOneCat){
		if (isOneCat){
			rdbtnCatYes.setSelected(true);
			rdbtnCatNo.setSelected(false);
			lblNumCat.setEnabled(false);
			cmbxNumCat.setEnabled(false);
			lblSiteRate.setEnabled(false);
			txtSiteRate.setEnabled(false);
			lblCatNum.setEnabled(false);
			cmbxRateCatNum.setEnabled(false);
			btnCatFile.setEnabled(false);
			txtCatFile.setEnabled(false);
		}
		else{
			rdbtnCatYes.setSelected(false);
			rdbtnCatNo.setSelected(true);
			lblNumCat.setEnabled(true);
			cmbxNumCat.setEnabled(true);
			lblSiteRate.setEnabled(true);
			txtSiteRate.setEnabled(true);
			lblCatNum.setEnabled(true);
			cmbxRateCatNum.setEnabled(true);
			btnCatFile.setEnabled(true);
			txtCatFile.setEnabled(true);
		}
	}

	
	protected void SitesWeightToggle(boolean isSites){
		if (isSites){
			rbtnPosWeightYes.setSelected(false);
			rbtnPosWeightNo.setSelected(true);
			btnWeightFile.setEnabled(false);
			txtWeightFile.setEnabled(false);
		}
		else{
			rbtnPosWeightYes.setSelected(true);
			rbtnPosWeightNo.setSelected(false);
			btnWeightFile.setEnabled(true);
			txtWeightFile.setEnabled(true);
		}
	}
	
	
	protected void EqualBFToggle(boolean isEmp){
		if (isEmp){
			rdbtnEqualBFYes.setSelected(true);
			rdbtnEqualBFNo.setSelected(false);
			lblBaseFreq.setEnabled(false);
			lblBaseFreqA.setEnabled(false);
			txtBaseFreqA.setEnabled(false);
			lblBaseFreqC.setEnabled(false);
			txtBaseFreqC.setEnabled(false);
			lblBaseFreqG.setEnabled(false);
			txtBaseFreqG.setEnabled(false);
			lblBaseFreqTU.setEnabled(false);
			txtBaseFreqTU.setEnabled(false);
		}
		else{
			rdbtnEqualBFYes.setSelected(false);
			rdbtnEqualBFNo.setSelected(true);
			lblBaseFreq.setEnabled(true);
			lblBaseFreqA.setEnabled(true);
			txtBaseFreqA.setEnabled(true);
			lblBaseFreqC.setEnabled(true);
			txtBaseFreqC.setEnabled(true);
			lblBaseFreqG.setEnabled(true);
			txtBaseFreqG.setEnabled(true);
			lblBaseFreqTU.setEnabled(true);
			txtBaseFreqTU.setEnabled(true);
		}
	}
	
	protected void MultToggle(boolean isMult){
		if(isMult){
			rbtnMultDataNo.setSelected(false);
			rbtnMultDataYes.setSelected(true);
			rbtnDataSets.setEnabled(true);
			rbtnWeights.setEnabled(true);
			lblMultData.setEnabled(true);
			lblHowManyData.setEnabled(true);
			txtNumSeqs.setEnabled(true);
			lblInputSeq.setEnabled(true);
			rdbtnInputSeqYes.setEnabled(true);
			rdbtnInputSeqNo.setEnabled(true);
		}
		else{
			rbtnMultDataNo.setSelected(true);
			rbtnMultDataYes.setSelected(false);
			rbtnDataSets.setEnabled(false);
			rbtnWeights.setEnabled(false);
			lblMultData.setEnabled(false);
			lblHowManyData.setEnabled(false);
			txtNumSeqs.setEnabled(false);
			lblInputSeq.setEnabled(false);
			rdbtnInputSeqYes.setEnabled(false);
			rdbtnInputSeqNo.setEnabled(false);
		}
	}
	
	protected void InputSeqToggle(boolean isSeq){
		if(isSeq){
			rdbtnInputSeqYes.setSelected(true);
			rdbtnInputSeqNo.setSelected(false);
		}
		else{
			rdbtnInputSeqYes.setSelected(false);
			rdbtnInputSeqNo.setSelected(true);
		}
	}
	
	protected void PrintDataToggle(boolean isData){
		if(isData){
			rbtnPrintDataYes.setSelected(true);
			rbtnPrintDataNo.setSelected(false);
		}
		else
		{
			rbtnPrintDataYes.setSelected(false);
			rbtnPrintDataNo.setSelected(true);
		}
	}
	
	protected void PrintIndToggle(boolean isInd){
		if(isInd){
			rbtnPrintIndYes.setSelected(true);
			rbtnPrintIndNo.setSelected(false);
		}
		else{
			rbtnPrintIndYes.setSelected(false);
			rbtnPrintIndNo.setSelected(true);
		}
	}
	
	
	protected void DataWeightToggle(boolean isData){
		if(isData){
			rbtnDataSets.setSelected(true);
			rbtnWeights.setSelected(false);
			//RandOrderToggle(true);
		}
		else{
			rbtnDataSets.setSelected(false);
			rbtnWeights.setSelected(true);
			//RandOrderToggle(false);
		}
	}

	protected void DisplaySiteRateValue(int catnum){
		if ((catnum-1) != lastsiteratecat){
			SiteRateValues[lastsiteratecat] = Double.parseDouble(txtSiteRate.getText());
			String siteratestr = Double.toString(SiteRateValues[catnum-1]);
			txtSiteRate.setText(siteratestr);
		}
		lastsiteratecat = catnum-1;
	}

	protected void DotDiffToggle(boolean isUseDot) {
		if (isUseDot) {
			rbtnDotDiffYes.setSelected(true);
			rbtnDotDiffNo.setSelected(false);
		} else {
			rbtnDotDiffYes.setSelected(false);
			rbtnDotDiffNo.setSelected(true);
		}
	}

	/**
	 * Create the application.
	 */
	public ProtDistUserInterface(String[] args) {
		phylipCall = false;
		if (args.length > 0)
		{
			if (args[0].contains("Phylip"))
			{
				phylipCall = true;
			}
		}
		initialize();
	}
	
	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		
		SiteRateValues = new double[9];
		SiteRateValues[0] = 1.0;
		SiteRateValues[1] = 1.0;
		SiteRateValues[2] = 1.0;
		SiteRateValues[3] = 1.0;
		SiteRateValues[4] = 1.0;
		SiteRateValues[5] = 1.0;
		SiteRateValues[6] = 1.0;
		SiteRateValues[7] = 1.0;
		SiteRateValues[8] = 1.0;
		
		lastsiteratecat = 0;
		
		filedir = System.getProperty("user.dir");

		frmProtDistControls = new JFrame();
		frmProtDistControls.setBackground(new Color(204, 255, 255));
		frmProtDistControls.setTitle("Protdist");
		frmProtDistControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmProtDistControls.setBounds(100, 100, 800, 810);
		frmProtDistControls.setPreferredSize(new Dimension(frmProtDistControls.getBounds().width, frmProtDistControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmProtDistControls.getPreferredSize());
		frmProtDistControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmProtDistControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow]", "[][][][]"));
		
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
		txtInputFile.setText("infile");
		txtInputFile.setBounds(168, 11, 467, 20);
		panel.add(txtInputFile, "cell 1 0 2 1,growx");
	
		btnWeightFile = new JButton("Weights File");
		btnWeightFile.setEnabled(false);
		btnWeightFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightFile);
			}
		});
		btnWeightFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnWeightFile, "cell 0 1,growx");
	
		txtWeightFile = new JTextField();
		txtWeightFile.setEnabled(false);
		txtWeightFile.setText("weightfile");
		txtWeightFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtWeightFile, "cell 1 1 2 1,growx");
		
		btnCatFile = new JButton("Categories File");
		btnCatFile.setEnabled(false);
		btnCatFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtCatFile);
			}
		});
		btnCatFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnCatFile, "cell 0 2,growx");
		
		txtCatFile = new JTextField();
		txtCatFile.setEnabled(false);
		txtCatFile.setText("catfile");
		txtCatFile.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(txtCatFile, "cell 1 2 2 1,growx");
		
		btnOutputFile = new JButton("Output File");
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		panel.add(btnOutputFile, "cell 0 3,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 3 2 1,growx");
	
		lblModel = new JLabel("Distance Model:");
		lblModel.setHorizontalAlignment(SwingConstants.RIGHT);
		lblModel.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblModel, "flowx,cell 0 4 2 1,alignx right");
		
		cmbxModel = new JComboBox();
		cmbxModel.setModel(new DefaultComboBoxModel(new String[] {"Jones-Taylor-Thornton", "Henikoff/Tillier PMB", "Dayhoff PAM", 
																  "mtRev", "Whelan and Goldman (WAG)", "mtMam", "Kimura protein", 
																  "Similarity table", "Categories"}));
		cmbxModel.setSelectedIndex(0);
		cmbxModel.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ModelToggle(cmbxModel.getSelectedIndex());
			}
		});
		cmbxModel.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxModel, "cell 2 4, growx");
				
		lblGamma = new JLabel("Gamma distribution of weights:");
		lblGamma.setHorizontalAlignment(SwingConstants.RIGHT);
		lblGamma.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGamma, "flowx,cell 0 5 2 1,alignx right");
		
		cmbxGamma = new JComboBox();
		cmbxGamma.setModel(new DefaultComboBoxModel(new String[] {"No", "Yes", "Gamma + invariant"}));
		cmbxGamma.setSelectedIndex(0);
		cmbxGamma.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				GammaToggle(cmbxGamma.getSelectedIndex());
			}
		});
		panel.add(cmbxGamma, "cell 2 5, growx");
		
		lblFracInvar = new JLabel("Fraction of sites invariant:");
		lblFracInvar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblFracInvar.setFont(new Font("Arial", Font.BOLD, 13));
		lblFracInvar.setEnabled(false);
		panel.add(lblFracInvar, "flowx,cell 0 6 2 1,alignx right");
	
		txtFracInvar = new JTextField();
		txtFracInvar.setText("0.0");
		txtFracInvar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtFracInvar.setEnabled(false);
		txtFracInvar.setColumns(6);
		panel.add(txtFracInvar, "cell 2 6");
		
		lblCoeffVar = new JLabel("Coefficient of variation:");
		lblCoeffVar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCoeffVar.setFont(new Font("Arial", Font.BOLD, 13));
		lblCoeffVar.setEnabled(false);
		panel.add(lblCoeffVar, "flowx,cell 0 7 2 1,alignx right");

		txtCoeffVar = new JTextField();
		txtCoeffVar.setEnabled(false);
		txtCoeffVar.setText("1.0");
		txtCoeffVar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCoeffVar.setColumns(6);
		panel.add(txtCoeffVar, "cell 2 7");
		
		lblCoeffVarNote = new JLabel("(for gamma dist = 1/\u221Aalpha)");
		lblCoeffVarNote.setEnabled(false);
		lblCoeffVarNote.setHorizontalAlignment(SwingConstants.LEFT);
		lblCoeffVarNote.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCoeffVarNote, "cell 2 7");

		lblCatSites = new JLabel("One category of sites:");
		lblCatSites.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCatSites.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCatSites, "flowx,cell 0 8 2 1,alignx right");
		
		rdbtnCatYes = new JRadioButton("Yes");
		rdbtnCatYes.setSelected(true);
		rdbtnCatYes.setBackground(new Color(204, 255, 255));
		rdbtnCatYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCatYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CatToggle(true);
			}
		});
		panel.add(rdbtnCatYes, "cell 2 8");
		
		rdbtnCatNo = new JRadioButton("No");
		rdbtnCatNo.setBackground(new Color(204, 255, 255));
		rdbtnCatNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnCatNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CatToggle(false);
			}
		});
		panel.add(rdbtnCatNo, "cell 2 8");
		
		lblNumCat = new JLabel("Number of site categories:");
		lblNumCat.setEnabled(false);
		lblNumCat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumCat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNumCat, "flowx,cell 0 9 2 1,alignx right");
		
		cmbxNumCat = new JComboBox();
		cmbxNumCat.setEnabled(false);
		cmbxNumCat.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxNumCat.setSelectedIndex(1);
		cmbxNumCat.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(cmbxNumCat, "cell 2 9");
		
		lblCatNum = new JLabel("Category:");
		lblCatNum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCatNum.setFont(new Font("Arial", Font.BOLD, 13));
		lblCatNum.setEnabled(false);
		panel.add(lblCatNum, "cell 2 9");
		
		cmbxRateCatNum = new JComboBox();
		cmbxRateCatNum.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cmbxRateCatNum.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxRateCatNum.setEnabled(false);
		cmbxRateCatNum.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplaySiteRateValue(Integer.parseInt(cmbxRateCatNum.getSelectedItem().toString()));
			}
		});
		panel.add(cmbxRateCatNum, "cell 2 9");
		
		lblSiteRate = new JLabel("Rate:");
		lblSiteRate.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSiteRate.setFont(new Font("Arial", Font.BOLD, 13));
		lblSiteRate.setEnabled(false);
		panel.add(lblSiteRate, "cell 2 9");
		
		txtSiteRate = new JTextField();
		txtSiteRate.setText("1.0");
		txtSiteRate.setHorizontalAlignment(SwingConstants.CENTER);
		txtSiteRate.setFont(new Font("Arial", Font.PLAIN, 13));
		txtSiteRate.setEnabled(false);
		txtSiteRate.setColumns(6);
		panel.add(txtSiteRate, "cell 2 9");
		
		lblPosWeight = new JLabel("Positions weighted:");
		lblPosWeight.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPosWeight.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPosWeight, "flowx,cell 0 10 2 1,alignx right");
		
		rbtnPosWeightYes = new JRadioButton("Yes");
		rbtnPosWeightYes.setBackground(new Color(204, 255, 255));
		rbtnPosWeightYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPosWeightYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SitesWeightToggle(false);
			}
		});
		rbtnPosWeightYes.setSelected(false);
		panel.add(rbtnPosWeightYes, "cell 2 10");
		
		rbtnPosWeightNo = new JRadioButton("No");
		rbtnPosWeightNo.setBackground(new Color(204, 255, 255));
		rbtnPosWeightNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPosWeightNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SitesWeightToggle(true);
			}
		});
		rbtnPosWeightNo.setSelected(true);
		panel.add(rbtnPosWeightNo, "cell 2 10");
		
		lblGenCode = new JLabel("Genetic code to use:");
		lblGenCode.setEnabled(false);
		lblGenCode.setHorizontalAlignment(SwingConstants.RIGHT);
		lblGenCode.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGenCode, "flowx,cell 0 11 2 1,alignx right");
		
		cmbxGenCode = new JComboBox();
		cmbxGenCode.setEnabled(false);
		cmbxGenCode.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxGenCode.setModel(new DefaultComboBoxModel(new String[] {"Universal", "Mitochondrial", "Vertebrate mitochondrial", "Fly mitochondrial", "Yeast mitochondrial"}));
		cmbxGenCode.setSelectedIndex(0);
		panel.add(cmbxGenCode, "cell 2 11");
		
		lblAACat = new JLabel("Amino Acid catigorization:");
		lblAACat.setEnabled(false);
		lblAACat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAACat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAACat, "flowx,cell 0 12 2 1,alignx right");
			
		lblAAInfo = new JLabel("[All have groups: (Glu Gln Asp Asn), (Lys Arg His), (Phe Tyr Trp)]");
		lblAAInfo.setEnabled(false);
		lblAAInfo.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAAInfo.setFont(new Font("Arial", Font.PLAIN, 13));
		panel.add(lblAAInfo, "cell 2 12");
		
		lblAddAA = new JLabel("Add:");
		lblAddAA.setEnabled(false);
		lblAddAA.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAddAA.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAddAA, "flowx,cell 0 13 2 1,alignx right");

		cmbxAACat = new JComboBox();
		cmbxAACat.setEnabled(false);
		cmbxAACat.setFont(new Font("Arial", Font.PLAIN, 13));
		cmbxAACat.setModel(new DefaultComboBoxModel(new String[] {"George/Hunt/Barker: (Cys), (Met Val Leu Ileu), (Gly  Ala  Ser  Thr Pro)", 
																  "Chemical: (Cys Met), (Val Leu Ileu Gly Ala Ser Thr), (Pro)", 
																  "Hall: (Cys), (Met Val Leu Ileu), (Gly Ala Ser Thr), (Pro)"}));
		cmbxAACat.setSelectedIndex(0);
		panel.add(cmbxAACat, "cell 2 13");
		
		lblProbChange = new JLabel("Probability of category change:");
		lblProbChange.setHorizontalAlignment(SwingConstants.RIGHT);
		lblProbChange.setFont(new Font("Arial", Font.BOLD, 13));
		lblProbChange.setEnabled(false);
		panel.add(lblProbChange, "flowx,cell 0 14 2 1,alignx right");
		
		txtProbChange = new JTextField();
		txtProbChange.setText("0.4570");
		txtProbChange.setFont(new Font("Arial", Font.PLAIN, 13));
		txtProbChange.setEnabled(false);
		txtProbChange.setColumns(6);
		panel.add(txtProbChange, "cell 2 14");
		
		lblProbChangeInfo = new JLabel("(1.0 = easy)");
		lblProbChangeInfo.setHorizontalAlignment(SwingConstants.RIGHT);
		lblProbChangeInfo.setFont(new Font("Arial", Font.PLAIN, 13));
		lblProbChangeInfo.setEnabled(false);
		panel.add(lblProbChangeInfo, "cell 2 14");
		
		lblTTRatio = new JLabel("Transition/transversion ratio:");
		lblTTRatio.setEnabled(false);
		lblTTRatio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblTTRatio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTTRatio, "flowx,cell 0 15 2 1,alignx right");
		
		txtTTRatio = new JTextField();
		txtTTRatio.setEnabled(false);
		txtTTRatio.setText("2.00");
		txtTTRatio.setFont(new Font("Arial", Font.PLAIN, 13));
		txtTTRatio.setColumns(6);
		panel.add(txtTTRatio, "cell 2 15");
		
		lblEqualBF = new JLabel("Use equal base frequencies:");
		lblEqualBF.setEnabled(false);
		lblEqualBF.setHorizontalAlignment(SwingConstants.RIGHT);
		lblEqualBF.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblEqualBF, "flowx,cell 0 16 2 1,alignx right");
		
		rdbtnEqualBFYes = new JRadioButton("Yes");
		rdbtnEqualBFYes.setEnabled(false);
		rdbtnEqualBFYes.setSelected(true);
		rdbtnEqualBFYes.setBackground(new Color(204, 255, 255));
		rdbtnEqualBFYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnEqualBFYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EqualBFToggle(true);
			}
		});
		panel.add(rdbtnEqualBFYes, "cell 2 16");
		
		rdbtnEqualBFNo = new JRadioButton("No");
		rdbtnEqualBFNo.setEnabled(false);
		rdbtnEqualBFNo.setBackground(new Color(204, 255, 255));
		rdbtnEqualBFNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnEqualBFNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				EqualBFToggle(false);
			}
		});
		panel.add(rdbtnEqualBFNo, "cell 2 16");
		
		lblBaseFreq = new JLabel("Base frequencies:");
		lblBaseFreq.setEnabled(false);
		lblBaseFreq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblBaseFreq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreq, "flowx,cell 0 17 2 1,alignx right");
		
		lblBaseFreqA = new JLabel("A");
		lblBaseFreqA.setEnabled(false);
		lblBaseFreqA.setBackground(new Color(153, 255, 255));
		lblBaseFreqA.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqA, "cell 2 17");
		
		txtBaseFreqA = new JTextField();
		txtBaseFreqA.setEnabled(false);
		txtBaseFreqA.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqA.setText("0.25");
		txtBaseFreqA.setColumns(5);
		panel.add(txtBaseFreqA, "cell 2 17");
		
		lblBaseFreqC = new JLabel("C");
		lblBaseFreqC.setEnabled(false);
		lblBaseFreqC.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqC, "cell 2 17");
		
		txtBaseFreqC = new JTextField();
		txtBaseFreqC.setEnabled(false);
		txtBaseFreqC.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqC.setText("0.25");
		txtBaseFreqC.setColumns(5);
		panel.add(txtBaseFreqC, "cell 2 17");
		
		lblBaseFreqG = new JLabel("G");
		lblBaseFreqG.setEnabled(false);
		lblBaseFreqG.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqG, "cell 2 17");
		
		txtBaseFreqG = new JTextField();
		txtBaseFreqG.setEnabled(false);
		txtBaseFreqG.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqG.setText("0.25");
		txtBaseFreqG.setColumns(5);
		panel.add(txtBaseFreqG, "cell 2 17");
		
		lblBaseFreqTU = new JLabel("T/U");
		lblBaseFreqTU.setEnabled(false);
		lblBaseFreqTU.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreqTU, "cell 2 17");
		
		txtBaseFreqTU = new JTextField();
		txtBaseFreqTU.setEnabled(false);
		txtBaseFreqTU.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqTU.setText("0.25");
		txtBaseFreqTU.setColumns(5);
		panel.add(txtBaseFreqTU, "cell 2 17");
		
		lblAnalyzeMultData = new JLabel("Analyze multiple data sets:");
		lblAnalyzeMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblAnalyzeMultData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblAnalyzeMultData, "flowx,cell 0 18 2 1,alignx right");
		
		rbtnMultDataYes = new JRadioButton("Yes");
		rbtnMultDataYes.setBackground(new Color(204, 255, 255));
		rbtnMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		rbtnMultDataYes.setSelected(false);
		panel.add(rbtnMultDataYes, "cell 2 18");
		
		rbtnMultDataNo = new JRadioButton("No");
		rbtnMultDataNo.setBackground(new Color(204, 255, 255));
		rbtnMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		rbtnMultDataNo.setSelected(true);
		panel.add(rbtnMultDataNo, "cell 2 18");
		
		lblMultData = new JLabel("Multiple data sets or multiple weights:");
		lblMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		lblMultData.setEnabled(false);
		panel.add(lblMultData, "flowx,cell 0 19 2 1,alignx right");
		
		rbtnDataSets = new JRadioButton("Data sets");
		rbtnDataSets.setBackground(new Color(204, 255, 255));
		rbtnDataSets.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(true);
			}
		});
		rbtnDataSets.setSelected(true);
		rbtnDataSets.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnDataSets.setEnabled(false);
		panel.add(rbtnDataSets, "cell 2 19");
		
		rbtnWeights = new JRadioButton("Weights");
		rbtnWeights.setBackground(new Color(204, 255, 255));
		rbtnWeights.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(false);
			}
		});
		
		rbtnWeights.setSelected(false);
		rbtnWeights.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnWeights.setEnabled(false);
		panel.add(rbtnWeights, "cell 2 19");

		lblHowManyData = new JLabel("Number:");
		lblHowManyData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHowManyData.setFont(new Font("Arial", Font.BOLD, 13));
		lblHowManyData.setEnabled(false);
		panel.add(lblHowManyData, "cell 2 19");
		
		txtNumSeqs = new JTextField();
		txtNumSeqs.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNumSeqs.setText("1");
		txtNumSeqs.setEnabled(false);
		txtNumSeqs.setColumns(6);
		panel.add(txtNumSeqs, "cell 2 19");
		
		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setEnabled(false);
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblInputSeq, "flowx,cell 0 20 2 1,alignx right");
		
		rdbtnInputSeqYes = new JRadioButton("Interleaved");
		rdbtnInputSeqYes.setEnabled(false);
		rdbtnInputSeqYes.setSelected(true);
		rdbtnInputSeqYes.setBackground(new Color(204, 255, 255));
		rdbtnInputSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		panel.add(rdbtnInputSeqYes, "cell 2 20");
		
		rdbtnInputSeqNo = new JRadioButton("Sequential");
		rdbtnInputSeqNo.setEnabled(false);
		rdbtnInputSeqNo.setBackground(new Color(204, 255, 255));
		rdbtnInputSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rdbtnInputSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		rdbtnInputSeqNo.setSelected(false);
		panel.add(rdbtnInputSeqNo, "cell 2 20");
		
		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintData, "flowx,cell 0 21 2 1,alignx right");
		
		rbtnPrintDataYes = new JRadioButton("Yes");
		rbtnPrintDataYes.setBackground(new Color(204, 255, 255));
		rbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		rbtnPrintDataYes.setSelected(false);
		panel.add(rbtnPrintDataYes, "cell 2 21");
		
		rbtnPrintDataNo = new JRadioButton("No");
		rbtnPrintDataNo.setBackground(new Color(204, 255, 255));
		rbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rbtnPrintDataNo.setSelected(true);
		panel.add(rbtnPrintDataNo, "cell 2 21");
		
		lblDotDiff = new JLabel("Use dot-difference display:");
		lblDotDiff.setFont(new Font("Arial", Font.BOLD, 13));
		lblDotDiff.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblDotDiff, "flowx,cell 0 22 2 1,alignx right");

		rbtnDotDiffYes = new JRadioButton("Yes");
		rbtnDotDiffYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnDotDiffYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDotDiffYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(true);
			}
		});
		rbtnDotDiffYes.setSelected(true);
		panel.add(rbtnDotDiffYes, "cell 2 22");

		rbtnDotDiffNo = new JRadioButton("No");
		rbtnDotDiffNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnDotDiffNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDotDiffNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(false);
			}
		});
		panel.add(rbtnDotDiffNo, "cell 2 22");
		
		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInd, "flowx,cell 0 23 2 1,alignx right");
		
		rbtnPrintIndYes = new JRadioButton("Yes");
		rbtnPrintIndYes.setBackground(new Color(204, 255, 255));
		rbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rbtnPrintIndYes.setSelected(true);
		panel.add(rbtnPrintIndYes, "cell 2 23");
		
		rbtnPrintIndNo = new JRadioButton("No");
		rbtnPrintIndNo.setBackground(new Color(204, 255, 255));
		rbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rbtnPrintIndNo.setSelected(false);
		panel.add(rbtnPrintIndNo, "cell 2 23");
		
		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			
			public void actionPerformed(ActionEvent e) {
				// catch last array data entries
				SiteRateValues[Integer.parseInt(cmbxRateCatNum.getSelectedItem().toString())-1] = Double.parseDouble(txtSiteRate.getText());

				inputvals = new ProtDistData();
				inputvals.infile = txtInputFile.getText();
				inputvals.wgtfile = txtWeightFile.getText();
				inputvals.catsfile = txtCatFile.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				
				if (cmbxModel.getSelectedIndex() == 0)
				{
					inputvals.Model = "JTT";
				}
				else if (cmbxModel.getSelectedIndex() == 1)
				{
					inputvals.Model = "HTPMB";
				}
				else if (cmbxModel.getSelectedIndex() == 2)
				{
					inputvals.Model = "DPAM";
				}
				else if (cmbxModel.getSelectedIndex() == 3)
				{
					inputvals.Model = "mtRev";
				}
				else if (cmbxModel.getSelectedIndex() == 4)
				{
					inputvals.Model = "WAG";
				}
				else if (cmbxModel.getSelectedIndex() == 5)
				{
					inputvals.Model = "mtMam";
				}
				else if (cmbxModel.getSelectedIndex() == 6)
				{
					inputvals.Model = "Kimura";
				}
				else if (cmbxModel.getSelectedIndex() == 7)
				{
					inputvals.Model = "Similarity";
				}
				else //(cmbxModel.getSelectedIndex() == 8)
				{
					inputvals.Model = "Cat";
				}
				
				if (cmbxGamma.getSelectedIndex() == 0)
				{
					inputvals.GammaDist = "No";
				}
				else if (cmbxGamma.getSelectedIndex() == 1)
				{
					inputvals.GammaDist = "Gamma";
				}
				else//  (cmbxGamma.getSelectedIndex() == 2)
				{
					inputvals.GammaDist = "Invariant";
				}
				inputvals.FracInvar = Double.parseDouble(txtFracInvar.getText());
				inputvals.CoeffVar = Double.parseDouble(txtCoeffVar.getText());
				
				inputvals.OneSite  = rdbtnCatYes.isSelected();
				inputvals.NumSites = Integer.parseInt(cmbxNumCat.getSelectedItem().toString());
				inputvals.SiteRate1 = SiteRateValues[0];
				inputvals.SiteRate2 = SiteRateValues[1];
				inputvals.SiteRate3 = SiteRateValues[2];
				inputvals.SiteRate4 = SiteRateValues[3];
				inputvals.SiteRate5 = SiteRateValues[4];
				inputvals.SiteRate6 = SiteRateValues[5];
				inputvals.SiteRate7 = SiteRateValues[6];
				inputvals.SiteRate8 = SiteRateValues[7];
				inputvals.SiteRate9 = SiteRateValues[8];
				
				inputvals.SitesWeighted = rbtnPosWeightYes.isSelected();

				if (cmbxGenCode.getSelectedIndex() == 0)
				{
					inputvals.GeneCode = "Universal";
				}
				else if (cmbxGenCode.getSelectedIndex() == 1)
				{
					inputvals.GeneCode = "Mitochondrial";
				}
				else if (cmbxGenCode.getSelectedIndex() == 2)
				{
					inputvals.GeneCode = "Vertebrate";
				}
				else if (cmbxGenCode.getSelectedIndex() == 3)
				{
					inputvals.GeneCode = "Fly";
				}
				else // if (cmbxGenCode.getSelectedIndex() == 4)
				{
					inputvals.GeneCode = "Yeast";
				}
				
				if (cmbxAACat.getSelectedIndex() == 0)
				{
					inputvals.AACat = "GHB";
				}
				else if (cmbxAACat.getSelectedIndex() == 1)
				{
					inputvals.AACat = "Chemical";
				}
				else //if (cmbxAACat.getSelectedIndex() == 2)
				{
					inputvals.AACat = "Hall";
				}
				inputvals.ProbCat = Double.parseDouble(txtProbChange.getText());
			
				inputvals.TTratio = Double.parseDouble(txtTTRatio.getText());
				inputvals.useEqualBF = rdbtnEqualBFYes.isSelected();
				inputvals.BaseFreqA = Double.parseDouble(txtBaseFreqA.getText());
				inputvals.BaseFreqC = Double.parseDouble(txtBaseFreqC.getText());
				inputvals.BaseFreqG = Double.parseDouble(txtBaseFreqG.getText());
				inputvals.BaseFreqTU = Double.parseDouble(txtBaseFreqTU.getText());

				inputvals.MultData = rbtnMultDataYes.isSelected();
				inputvals.MultDSet = rbtnDataSets.isSelected();
				inputvals.NumSets = Integer.parseInt(txtNumSeqs.getText());
				inputvals.InputSeq = rdbtnInputSeqYes.isSelected();
				inputvals.PrintData = rbtnPrintDataYes.isSelected();
				inputvals.DotDiff = rbtnDotDiffYes.isSelected();
				inputvals.PrintInd = rbtnPrintIndYes.isSelected();
					
				
				btnExecute.setEnabled(false);	
				String title = "Protdist Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread protDistThread = new Thread() {
						public void run() {
							runProtDistThreads();
						}
			  	    };
			  	    protDistThread.start();
				}
				btnExecute.setEnabled(true);
	
				}
			
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 2 24,alignx center");
		
		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				if(phylipCall)
				{
					frmProtDistControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 24,alignx center");
	

	}
	public boolean checkInputVals(){
		
		// check files
		TestFileNames test = new TestFileNames();
		
		if (!test.DuplicateFileNames(inputvals.infile, "Input", inputvals.outfile, "Output"))
		{			
			return false;		
		}

		if (!test.FileAvailable(inputvals.infile, "Input"))
		{
			return false;
		}

		String opt = test.FileAlreadyExists(inputvals.outfile, "Outfile");
		if (opt == "q")
		{
			return false;
		}
		else
		{
			inputvals.outfileopt = opt;
		}
		
		if (inputvals.SitesWeighted){
			if (!test.FileAvailable(inputvals.wgtfile, "Weights"))
			{
				return false;
			}
		}
		
		if (!inputvals.OneSite){
			if (!test.FileAvailable(inputvals.catsfile, "Categories"))
			{
				return false;
			}
		}

		// check data
		if (inputvals.BaseFreqA < 0) {
			String msg1 = "Input value: Base Frequency for A cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.BaseFreqC < 0) {
			String msg1 = "Input value: Base Frequency for C cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.BaseFreqG < 0) {
			String msg1 = "Input value: Base Frequency for G cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.BaseFreqTU < 0) {
			String msg1 = "Input value: Base Frequency for T/U cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.TTratio < 0) {
			String msg1 = "Input value: Transition / Transversion Ratio cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.CoeffVar <= 0.0) {
			String msg1 = "Input value: Coefficient of Variation must be greater than 0.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.ProbCat < 0.0) {
			String msg1 = "Input value: Probability of category change cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;			
		}
		
		if (inputvals.ProbCat > 1.0) {
			String msg1 = "Input value: Probability of category change cannot be greater than 1.0";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;			
		}
		
		if (inputvals.FracInvar < 0.0) {
			String msg1 = "Input value: Fraction of sites invariant cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.FracInvar > 1.0) {
			String msg1 = "Input value: Fraction of sites invariant cannot greater than 1.0.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		// autoscale frequencies
		double sum = inputvals.BaseFreqA + inputvals.BaseFreqC + inputvals.BaseFreqG + inputvals.BaseFreqTU;
		inputvals.BaseFreqA = inputvals.BaseFreqA / sum;
		inputvals.BaseFreqC = inputvals.BaseFreqC / sum;
		inputvals.BaseFreqG = inputvals.BaseFreqG / sum;
		inputvals.BaseFreqTU = inputvals.BaseFreqTU / sum;
		return true;
	}
	
	protected void runProtDistThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("protdist", ProtDist.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("ProtDist");
    		return;
    	}
		try 
		{
	  	    Thread protDistRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			ProtDist protdist = (ProtDist) Native.loadLibrary("protdist", ProtDist.class);
		  			protdist.protdist(		
		  					inputvals.infile,
		  					inputvals.wgtfile,
		  					inputvals.catsfile,
		  					inputvals.outfile,
		  					inputvals.outfileopt,
		  					inputvals.Model,
		  					inputvals.GammaDist,
		  					inputvals.FracInvar,
		  					inputvals.CoeffVar,
		  					inputvals.OneSite,
		  					inputvals.NumSites,
		  					inputvals.SiteRate1,  
		  					inputvals.SiteRate2,  
		  					inputvals.SiteRate3,  
		  					inputvals.SiteRate4,  
		  					inputvals.SiteRate5,  
		  					inputvals.SiteRate6,  
		  					inputvals.SiteRate7,  
		  					inputvals.SiteRate8,  
		  					inputvals.SiteRate9,  
		  					inputvals.SitesWeighted,
		  					inputvals.GeneCode,
		  					inputvals.AACat,
		  					inputvals.ProbCat,
		  					inputvals.TTratio,
		  					inputvals.useEqualBF,
		  					inputvals.BaseFreqA,
		  					inputvals.BaseFreqC,
		  					inputvals.BaseFreqG,
		  					inputvals.BaseFreqTU,
		  					inputvals.MultData,
		  					inputvals.MultDSet,
		  					inputvals.NumSets,
		  					inputvals.InputSeq,
		  					inputvals.PrintData,
		  					inputvals.DotDiff,
		  					inputvals.PrintInd);
			  	    };
	  	    };
	  	    protDistRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (protDistRunThread.isAlive());
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
