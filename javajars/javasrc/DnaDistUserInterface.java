package phylip;

import java.awt.EventQueue;
import javax.swing.*;

import com.sun.jna.Library;
import com.sun.jna.Native;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import java.awt.Font;
import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class DnaDistUserInterface {
	
	   public interface DnaDist extends Library {
	        public void dnadist(
	        		String infile,
	        		String outfile,
	        		String outfileopt,
	        		String catfile,
	        		String weightfile,
	        		String distanceoptions,
	        		String gammadistributed,
	        		double coefvar,
	        		double pctinvar,
	        		boolean usettratio,
	        		double ttratio,
	        		boolean usesubstrates,
	        		int numcats,
	        		double rate1,
	        		double rate2,
	        		double rate3,
	        		double rate4,
	        		double rate5,
	        		double rate6,
	        		double rate7,
	        		double rate8,
	        		double rate9,
	        		boolean useweights,
	        		boolean useEmpBF,
	        		double basefreqA,
	        		double basefreqC,
	        		double basefreqG,
	        		double basefreqTU,
	        		String distmatrix,
	        		boolean usemultdataset,
	        		boolean useData,
	        		int numdatasets,
	        		boolean doseqint,
	        		boolean doprintdata,
	        		boolean dodotdif,
	        		boolean doprintind);
	    }

	public class DnaDistData {
		String infile;
		String outfile;
		String outfileopt;
		String catfile;
		String weightfile;
		String distanceoptions;
		String gammadistributed;
		double coefvar;
		double pctinvar;
		boolean usettratio;
		double ttratio;
		boolean usesubstrates;
		int numcats;
		boolean useweights;
		// these are explicitly named because JNA doesn't pass arrays gracefully
		double rate1;  
		double rate2;  
		double rate3;  
		double rate4;  
		double rate5;  
		double rate6;  
		double rate7;  
		double rate8;  
		double rate9;  
		boolean useEmpBF;
		double basefreqA;
		double basefreqC;
		double basefreqG;
		double basefreqTU;
		String distmatrix;
		boolean usemultdataset;
		boolean useData;
		int dataset;
		boolean doseqint;
		boolean doprintdata;
		boolean dotdiff;
		boolean doprintind;
	}

	private DnaDistData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	public enum Gammatype{YES, NO, GAMMA}
	private boolean ExplicitWgts;
	private boolean phylipCall;

	private JFrame frmDnaDistControls;
	private JLabel lblDistanceOptions;
	private JLabel lblGammaDistributedRates;
	private JRadioButton btnGammaNo;
	private JRadioButton btnGammaYes;
	private JRadioButton btnGammaInvariant;
	private JButton btnInFile;
	private JTextField txtInFile;
	private JTextField txtTTRatio;
	private JLabel lblTTRatio;
	private JLabel lblSubRates;
	private JRadioButton btnSubsYes;
	private JRadioButton btnSubsNo;
	private JLabel lblNumCats;
	//private JComboBox<String> cbxNumCats;
	private JComboBox cbxNumCats;
	private JLabel lblUseWeightsFor;
	private JRadioButton btnWeightYes;
	private JRadioButton btnWeightNo;
	private JLabel lblUseEmpiricalBase;
	private JRadioButton btnEmpYes;
	private JRadioButton btnEmpNo;
	private JLabel lblBaseFreq;
	private JLabel lblBaseFreqA;
	private JLabel lblBaseFreqC;
	private JLabel lblBaseFreqG;
	private JLabel lblBaseFreqTU;
	private JTextField txtBaseFreqA;
	private JTextField txtBaseFreqC;
	private JTextField txtBaseFreqG;
	private JTextField txtBaseFreqTU;
	private JLabel lblDistMat;
	//private JComboBox<Object> cbxDistMat;
	private JComboBox cbxDistMat;
	private JLabel lblMultData;
	private JRadioButton btnMultYes;
	private JRadioButton btnMultNo;
	private JLabel lblKindData;
	private JRadioButton btnMultDataSet;
	private JRadioButton btnMultWeight;
	private JLabel lblHowManyData;
	private JTextField txtDataWeight;
	private JLabel lblInputSeq;
	private JRadioButton btnInputSeqYes;
	private JRadioButton btnInputSeqNo;
	private JLabel lblPrintData;
	private JRadioButton btnPrintDataYes;
	private JRadioButton btnPrintDataNo;
	private JLabel lblDoProgress;
	private JRadioButton btnDoProgressYes;
	private JRadioButton btnDoProgressNo;
	private JButton btnExecute;
	private JButton btnQuit;
	private JButton btnOutFile;
	private JTextField txtOutFile;
	//private JComboBox<String> cbxDistOpts;
	private JComboBox cbxDistOpts;
	private JButton btnCatsFile;
	private JTextField txtCatsFile;
	private JButton btnWeightsFile;
	private JTextField txtWeightsFile;
	public String GammaValue;
	private JTextField txtCoefVar;
	private JLabel lblFractionOfSites;
	private JTextField txtPcntInvar;
	private JLabel lblCoefVar;
	private JLabel lblRate;
	private JTextField txtRate;
	private JLabel lblCatNum;
	private JComboBox cbxRateCatNum;
	private double[] ratevalues;
	private int lastratecat;
	private JLabel lblRatio;
	private JRadioButton rbtnTTYes;
	private JRadioButton rbtnTTNo;
	private JLabel lblDotDiff;
	private JRadioButton btnDotDiffYes;
	private JRadioButton btnDotDiffNo;
	
	private JScrollPane scrollPane;
	private JPanel panel;

	/**
	 * Launch the application.
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DnaDistUserInterface window = new DnaDistUserInterface(args);
					window.frmDnaDistControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	protected void ChooseFile(JTextField file) {		
		//Construct a new file choose whose default path is the path to this executable, which 
		//is returned by System.getProperty("user.dir")
		
		JFileChooser fileChooser = new JFileChooser(filedir);

		int option = fileChooser.showOpenDialog(frmDnaDistControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}

	/**
	 * Create the application.
	 */
	public DnaDistUserInterface(String[] args) {
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
	
	protected void GammaToggle(Gammatype value) {
		if (value == Gammatype.YES){
			btnGammaYes.setSelected(true);
			btnGammaNo.setSelected(false);
			btnGammaInvariant.setSelected(false);
			GammaValue = "Yes";
			lblCoefVar.setEnabled(true);
			txtCoefVar.setEnabled(true);
			lblFractionOfSites.setEnabled(false);
			txtPcntInvar.setEnabled(false);
			}
		else if (value == Gammatype.NO) {
			btnGammaYes.setSelected(false);
			btnGammaNo.setSelected(true);
			btnGammaInvariant.setSelected(false);	
			GammaValue = "No";
			lblCoefVar.setEnabled(false);
			txtCoefVar.setEnabled(false);
			lblFractionOfSites.setEnabled(false);
			txtPcntInvar.setEnabled(false);
			}
		else if (value == Gammatype.GAMMA) {
			btnGammaInvariant.setSelected(true);
			btnGammaYes.setSelected(false);
			btnGammaNo.setSelected(false);
			GammaValue = "Gamma";
			lblCoefVar.setEnabled(true);
			txtCoefVar.setEnabled(true);
			lblFractionOfSites.setEnabled(true);
			txtPcntInvar.setEnabled(true);
		}
	}

	protected void CategoryToggle(boolean iscat) {
		if (iscat){
			btnSubsYes.setSelected(true);
			btnSubsNo.setSelected(false);
		}
		else {
			btnSubsNo.setSelected(true);
			btnSubsYes.setSelected(false);
		}
	}
	
	protected void WeightToggle(boolean isWeight) {
		if (isWeight){
			btnWeightYes.setSelected(true);
			btnWeightNo.setSelected(false);
			btnWeightsFile.setEnabled(true);
			txtWeightsFile.setEnabled(true);
			ExplicitWgts = true;
		}
		else {
			btnWeightNo.setSelected(true);
			btnWeightYes.setSelected(false);
			btnWeightsFile.setEnabled(false);
			txtWeightsFile.setEnabled(false);
			ExplicitWgts = false;
		}
	}
	
	protected void EmpToggle(boolean isEmp) {
		if (isEmp){
			btnEmpYes.setSelected(true);
			btnEmpNo.setSelected(false);
			lblBaseFreq.setEnabled(false);
			lblBaseFreqA.setEnabled(false);
			lblBaseFreqC.setEnabled(false);
			lblBaseFreqG.setEnabled(false);
			lblBaseFreqTU.setEnabled(false);
			txtBaseFreqA.setEnabled(false);
			txtBaseFreqC.setEnabled(false);
			txtBaseFreqG.setEnabled(false);
			txtBaseFreqTU.setEnabled(false);
		}
		else {
			btnEmpNo.setSelected(true);
			btnEmpYes.setSelected(false);
			lblBaseFreq.setEnabled(true);
			lblBaseFreqA.setEnabled(true);
			lblBaseFreqC.setEnabled(true);
			lblBaseFreqG.setEnabled(true);
			lblBaseFreqTU.setEnabled(true);
			txtBaseFreqA.setEnabled(true);
			txtBaseFreqC.setEnabled(true);
			txtBaseFreqG.setEnabled(true);
			txtBaseFreqTU.setEnabled(true);		
		}
	}
	
	protected void MultToggle(boolean isdata) {
		if (isdata){
			btnMultYes.setSelected(true);
			btnMultNo.setSelected(false);
			lblKindData.setEnabled(true);
			btnMultDataSet.setEnabled(true);
			btnMultWeight.setEnabled(true);
			lblHowManyData.setEnabled(true);
			txtDataWeight.setEnabled(true);
			lblInputSeq.setEnabled(true);
			btnInputSeqYes.setEnabled(true);
			btnInputSeqNo.setEnabled(true);
			if (btnMultWeight.isSelected())
			{
				if(!ExplicitWgts)
				{
					btnWeightsFile.setEnabled(true);
					txtWeightsFile.setEnabled(true);	
					btnWeightYes.setSelected(true);
					btnWeightNo.setSelected(false);
				}
			}
		}
		else {
			btnMultNo.setSelected(true);
			btnMultYes.setSelected(false);
			lblKindData.setEnabled(false);
			btnMultDataSet.setEnabled(false);
			btnMultWeight.setEnabled(false);
			lblHowManyData.setEnabled(false);
			txtDataWeight.setEnabled(false);
			lblInputSeq.setEnabled(false);
			btnInputSeqYes.setEnabled(false);
			btnInputSeqNo.setEnabled(false);
			if(ExplicitWgts)
			{
				btnWeightsFile.setEnabled(true);
				txtWeightsFile.setEnabled(true);	
				btnWeightYes.setSelected(true);
				btnWeightNo.setSelected(false);
			}
			else
			{
				btnWeightsFile.setEnabled(false);
				txtWeightsFile.setEnabled(false);	
				btnWeightYes.setSelected(false);
				btnWeightNo.setSelected(true);				
			}
		}
	}
	
	protected void DataWeightToggle(boolean dataweight) {
		if (dataweight){
			btnMultDataSet.setSelected(true);
			btnMultWeight.setSelected(false);
			if(!ExplicitWgts)
			{
				btnWeightsFile.setEnabled(false);
				txtWeightsFile.setEnabled(false);
				btnWeightYes.setSelected(false);
				btnWeightNo.setSelected(true);
			}
	}
		else {
			btnMultWeight.setSelected(true);
			btnMultDataSet.setSelected(false);
			btnWeightsFile.setEnabled(true);
			txtWeightsFile.setEnabled(true);
			btnWeightYes.setSelected(true);
			btnWeightNo.setSelected(false);
		}
	}
	
	protected void SeqToggle(boolean isSeq) {
		if (isSeq){
			btnInputSeqYes.setSelected(true);
			btnInputSeqNo.setSelected(false);
		}
		else {
			btnInputSeqNo.setSelected(true);
			btnInputSeqYes.setSelected(false);
		}
	}
	
	protected void doprintdataToggle(boolean isPrint1) {
		if (isPrint1){
			btnPrintDataYes.setSelected(true);
			btnPrintDataNo.setSelected(false);
		}
		else {
			btnPrintDataNo.setSelected(true);
			btnPrintDataYes.setSelected(false);
		}
	}
	
	protected void doprintindToggle(boolean isPrint2) {
		if (isPrint2){
			btnDoProgressYes.setSelected(true);
			btnDoProgressNo.setSelected(false);
		}
		else {
			btnDoProgressNo.setSelected(true);
			btnDoProgressYes.setSelected(false);
		}
	}
	
	protected void DistToggle() {
		if ((cbxDistOpts.getSelectedItem().toString()).contains("F84")) {
			lblUseEmpiricalBase.setEnabled(true);
			btnEmpYes.setEnabled(true);
			btnEmpNo.setEnabled(true);
			if (btnEmpYes.isSelected())
			{
				lblBaseFreq.setEnabled(true);
				lblBaseFreqA.setEnabled(true);
				lblBaseFreqC.setEnabled(true);
				lblBaseFreqG.setEnabled(true);
				lblBaseFreqTU.setEnabled(true);
				txtBaseFreqA.setEnabled(true);
				txtBaseFreqC.setEnabled(true);
				txtBaseFreqG.setEnabled(true);
				txtBaseFreqTU.setEnabled(true);
			}
			else
			{
				lblBaseFreq.setEnabled(false);
				lblBaseFreqA.setEnabled(false);
				lblBaseFreqC.setEnabled(false);
				lblBaseFreqG.setEnabled(false);
				lblBaseFreqTU.setEnabled(false);
				txtBaseFreqA.setEnabled(false);
				txtBaseFreqC.setEnabled(false);
				txtBaseFreqG.setEnabled(false);
				txtBaseFreqTU.setEnabled(false);
				
			}
			lblTTRatio.setEnabled(true);
			txtTTRatio.setEnabled(true);
			lblGammaDistributedRates.setEnabled(true);
			btnGammaYes.setEnabled(true);
			btnGammaNo.setEnabled(true);
			btnGammaInvariant.setEnabled(true);
			lblSubRates.setEnabled(true);
			btnSubsYes.setEnabled(true);
			btnSubsNo.setEnabled(true);
		}
		else if ((cbxDistOpts.getSelectedItem().toString()).contains("Kimura")) {
			lblUseEmpiricalBase.setEnabled(false);
			btnEmpYes.setEnabled(false);
			btnEmpNo.setEnabled(false);
			lblBaseFreq.setEnabled(false);
			lblBaseFreqA.setEnabled(false);
			lblBaseFreqC.setEnabled(false);
			lblBaseFreqG.setEnabled(false);
			lblBaseFreqTU.setEnabled(false);
			txtBaseFreqA.setEnabled(false);
			txtBaseFreqC.setEnabled(false);
			txtBaseFreqG.setEnabled(false);
			txtBaseFreqTU.setEnabled(false);
			lblTTRatio.setEnabled(true);
			txtTTRatio.setEnabled(true);
			lblGammaDistributedRates.setEnabled(true);
			btnGammaYes.setEnabled(true);
			btnGammaNo.setEnabled(true);
			btnGammaInvariant.setEnabled(true);
			lblSubRates.setEnabled(true);
			btnSubsYes.setEnabled(true);
			btnSubsNo.setEnabled(true);
		}
		else if ((cbxDistOpts.getSelectedItem().toString()).contains("Jukes-Cantor")) {
			lblUseEmpiricalBase.setEnabled(false);
			btnEmpYes.setEnabled(false);
			btnEmpNo.setEnabled(false);
			lblBaseFreq.setEnabled(false);
			lblBaseFreqA.setEnabled(false);
			lblBaseFreqC.setEnabled(false);
			lblBaseFreqG.setEnabled(false);
			lblBaseFreqTU.setEnabled(false);
			txtBaseFreqA.setEnabled(false);
			txtBaseFreqC.setEnabled(false);
			txtBaseFreqG.setEnabled(false);
			txtBaseFreqTU.setEnabled(false);
			lblTTRatio.setEnabled(false);
			txtTTRatio.setEnabled(false);
			lblGammaDistributedRates.setEnabled(true);
			btnGammaYes.setEnabled(true);
			btnGammaNo.setEnabled(true);
			btnGammaInvariant.setEnabled(true);
			lblSubRates.setEnabled(true);
			btnSubsYes.setEnabled(true);
			btnSubsNo.setEnabled(true);
		}
		else if ((cbxDistOpts.getSelectedItem().toString()).contains("LogDet")) {
			lblUseEmpiricalBase.setEnabled(false);
			btnEmpYes.setEnabled(false);
			btnEmpNo.setEnabled(false);
			lblBaseFreq.setEnabled(false);
			lblBaseFreqA.setEnabled(false);
			lblBaseFreqC.setEnabled(false);
			lblBaseFreqG.setEnabled(false);
			lblBaseFreqTU.setEnabled(false);
			txtBaseFreqA.setEnabled(false);
			txtBaseFreqC.setEnabled(false);
			txtBaseFreqG.setEnabled(false);
			txtBaseFreqTU.setEnabled(false);
			lblTTRatio.setEnabled(false);
			txtTTRatio.setEnabled(false);
			lblGammaDistributedRates.setEnabled(false);
			btnGammaYes.setEnabled(false);
			btnGammaNo.setEnabled(false);
			btnGammaInvariant.setEnabled(false);
			lblSubRates.setEnabled(false);
			btnSubsYes.setEnabled(false);
			btnSubsNo.setEnabled(false);
			lblNumCats.setEnabled(false);
			cbxNumCats.setEnabled(false);
		}
		else {
			lblUseEmpiricalBase.setEnabled(false);
			btnEmpYes.setEnabled(false);
			btnEmpNo.setEnabled(false);
			lblBaseFreq.setEnabled(false);
			lblBaseFreqA.setEnabled(false);
			lblBaseFreqC.setEnabled(false);
			lblBaseFreqG.setEnabled(false);
			lblBaseFreqTU.setEnabled(false);
			txtBaseFreqA.setEnabled(false);
			txtBaseFreqC.setEnabled(false);
			txtBaseFreqG.setEnabled(false);
			txtBaseFreqTU.setEnabled(false);
			lblTTRatio.setEnabled(false);
			txtTTRatio.setEnabled(false);
			lblGammaDistributedRates.setEnabled(false);
			btnGammaYes.setEnabled(false);
			btnGammaNo.setEnabled(false);
			btnGammaInvariant.setEnabled(false);
			lblSubRates.setEnabled(false);
			btnSubsYes.setEnabled(false);
			btnSubsNo.setEnabled(false);
			lblNumCats.setEnabled(false);
			cbxNumCats.setEnabled(false);
		}
	}
	
	protected void CatsFileToggle(boolean catTog) {
		if (catTog) {
			btnSubsYes.setSelected(true);
			btnSubsNo.setSelected(false);
			btnCatsFile.setEnabled(true);
			txtCatsFile.setEnabled(true);
			lblNumCats.setEnabled(true);
			cbxNumCats.setEnabled(true);
			lblRate.setEnabled(true);
			txtRate.setEnabled(true);
			lblCatNum.setEnabled(true);
			cbxRateCatNum.setEnabled(true);
		}
		else {
			btnSubsNo.setSelected(true);
			btnSubsYes.setSelected(false);
			btnCatsFile.setEnabled(false);
			txtCatsFile.setEnabled(false);
			lblNumCats.setEnabled(false);
			cbxNumCats.setEnabled(false);
			lblRate.setEnabled(false);
			txtRate.setEnabled(false);
			lblCatNum.setEnabled(false);
			cbxRateCatNum.setEnabled(false);
			}
		}
		
	protected void CategoryFileToggle(){
		if ((cbxNumCats.getSelectedItem().toString()).contains("1")) {
			btnCatsFile.setEnabled(false);
			txtCatsFile.setEnabled(false);
		}
		else {
			btnCatsFile.setEnabled(true);
			txtCatsFile.setEnabled(true);
		}
	}

	protected void DisplayRateValue(int catnum){
		if ((catnum-1) != lastratecat){
			ratevalues[lastratecat] = Double.parseDouble(txtRate.getText());
			String ratestr = Double.toString(ratevalues[catnum-1]);
			txtRate.setText(ratestr);
		}
		lastratecat = catnum-1;
	}
	
	protected void 	EnableTTValue(boolean enableTTvalue){
		if(enableTTvalue){
			rbtnTTYes.setSelected(true);
			rbtnTTNo.setSelected(false);
			lblRatio.setEnabled(true);
			txtTTRatio.setEnabled(true);			
		}
		else{
			rbtnTTNo.setSelected(true);
			rbtnTTYes.setSelected(false);
			lblRatio.setEnabled(false);
			txtTTRatio.setEnabled(false);			
		}
	}

	protected void DotDiffToggle(boolean isDot) {
		if (isDot) {
			btnDotDiffYes.setSelected(true);
			btnDotDiffNo.setSelected(false);
		} else {
			btnDotDiffYes.setSelected(false);
			btnDotDiffNo.setSelected(true);
		}
	}

	/**
	 * Initialize the contents of the frame.
	 */
	
	
	private void initialize() {
		filedir = System.getProperty("user.dir");
		ExplicitWgts = false;

		ratevalues = new double[9];
		ratevalues[0] = 1.0;
		ratevalues[1] = 1.0;
		ratevalues[2] = 1.0;
		ratevalues[3] = 1.0;
		ratevalues[4] = 1.0;
		ratevalues[5] = 1.0;
		ratevalues[6] = 1.0;
		ratevalues[7] = 1.0;
		ratevalues[8] = 1.0;
		
		lastratecat = 0;
		
		GammaValue = "No";
	
		frmDnaDistControls = new JFrame();
		frmDnaDistControls.setBackground(new Color(204, 255, 255));
		frmDnaDistControls.setTitle("Dnadist");
		frmDnaDistControls.setFont(new Font("Arial", Font.BOLD, 13));
		frmDnaDistControls.setBounds(100, 100, 700, 790);
		frmDnaDistControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmDnaDistControls.setPreferredSize(new Dimension(frmDnaDistControls.getBounds().width, frmDnaDistControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmDnaDistControls.getPreferredSize());
		frmDnaDistControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmDnaDistControls.getPreferredSize());
		scrollPane.setViewportView(panel);
		panel.setLayout(new MigLayout("", "[pref!,grow][30.00,grow][pref!,grow]", "[][][][][]"));
		
		btnInFile = new JButton("Input File");
		btnInFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnInFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtInFile);
			}
		});
		panel.add(btnInFile, "cell 0 0,growx");
		
		txtInFile = new JTextField();
		txtInFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtInFile.setText("infile");
		panel.add(txtInFile, "cell 1 0 2 1,growx");
	
		btnCatsFile = new JButton("Categories File");
		btnCatsFile.setEnabled(false);
		btnCatsFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtCatsFile);;				
				}
			});
		btnCatsFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnCatsFile, "cell 0 1,growx");
		
		txtCatsFile = new JTextField();
		txtCatsFile.setEnabled(false);
		txtCatsFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCatsFile.setText("catfile");
		panel.add(txtCatsFile, "cell 1 1 2 1,growx");
		
		btnWeightsFile = new JButton("Weights File");
		btnWeightsFile.setEnabled(false);
		btnWeightsFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtWeightsFile);;				
				}
			});
		btnWeightsFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnWeightsFile, "cell 0 2,growx");
		
		txtWeightsFile = new JTextField();
		txtWeightsFile.setEnabled(false);
		txtWeightsFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtWeightsFile.setText("weightfile");
		panel.add(txtWeightsFile, "cell 1 2 2 1,growx");
		
		btnOutFile = new JButton("Output File");
		btnOutFile.setFont(new Font("Arial", Font.BOLD, 13));
		btnOutFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ChooseFile(txtOutFile);
			}
		});
		panel.add(btnOutFile, "cell 0 3,growx");
		
		txtOutFile = new JTextField();
		txtOutFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutFile.setText("outfile");
		panel.add(txtOutFile, "cell 1 3 2 1,growx");
		
		lblDistanceOptions = new JLabel("Distance Options:");
		lblDistanceOptions.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDistanceOptions.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDistanceOptions, "flowx,cell 0 4 2 1,alignx right");

		cbxDistOpts = new JComboBox();
		cbxDistOpts.setModel(new DefaultComboBoxModel(new String[] {"F84", "Kimura", "Jukes-Cantor", "LogDet", "Similarity Table"}));
		cbxDistOpts.setFont(new Font("Arial", Font.PLAIN, 13));
		cbxDistOpts.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DistToggle();
			}
		});
		panel.add(cbxDistOpts, "cell 2 4,growx");
	
		lblGammaDistributedRates = new JLabel("Gamma distributed rates across sites:");
		lblGammaDistributedRates.setHorizontalAlignment(SwingConstants.RIGHT);
		lblGammaDistributedRates.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGammaDistributedRates, "flowx,cell 0 5 2 1,alignx right");
		
		btnGammaYes = new JRadioButton("Yes");
		btnGammaYes.setFont(new Font("Arial", Font.BOLD, 13));
		btnGammaYes.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent arg0) {
			GammaToggle(Gammatype.YES);				
			}
		});
		btnGammaYes.setSelected(false);
		panel.add(btnGammaYes, "cell 2 5");
		
		btnGammaNo = new JRadioButton("No");
		btnGammaNo.setFont(new Font("Arial", Font.BOLD, 13));
		btnGammaNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				GammaToggle(Gammatype.NO);				
			}
		});
		btnGammaNo.setSelected(true);
		panel.add(btnGammaNo, "cell 2 5");
		
		btnGammaInvariant = new JRadioButton("Gamma+Invariant");
		btnGammaInvariant.setFont(new Font("Arial", Font.BOLD, 13));
		btnGammaInvariant.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GammaToggle(Gammatype.GAMMA);				
			}
		});
		btnGammaInvariant.setSelected(false);
		panel.add(btnGammaInvariant, "cell 2 5");
		
		lblCoefVar = new JLabel("Coefficient of Variation:");
		lblCoefVar.setEnabled(false);
		lblCoefVar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCoefVar.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCoefVar, "flowx,cell 0 6 2 1,alignx right");
		
		txtCoefVar = new JTextField();
		txtCoefVar.setEnabled(false);
		txtCoefVar.setText("1.0");
		txtCoefVar.setHorizontalAlignment(SwingConstants.CENTER);
		txtCoefVar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCoefVar.setColumns(6);
		panel.add(txtCoefVar, "cell 2 6");
		
		lblFractionOfSites = new JLabel("Fraction of Sites Invariant:");
		lblFractionOfSites.setEnabled(false);
		lblFractionOfSites.setHorizontalAlignment(SwingConstants.RIGHT);
		lblFractionOfSites.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblFractionOfSites, "flowx,cell 0 7 2 1,alignx right");
		
		txtPcntInvar = new JTextField();
		txtPcntInvar.setEnabled(false);
		txtPcntInvar.setText("0.0");
		txtPcntInvar.setHorizontalAlignment(SwingConstants.CENTER);
		txtPcntInvar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtPcntInvar.setColumns(6);
		panel.add(txtPcntInvar, "cell 2 7");
		
		lblTTRatio = new JLabel("Use Transition/Transversion ratio:");
		lblTTRatio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblTTRatio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTTRatio, "flowx,cell 0 8 2 1,alignx right");

		rbtnTTYes = new JRadioButton("Yes");
		rbtnTTYes.setSelected(false);
		rbtnTTYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnTTYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				EnableTTValue(true);
			}
		});
		panel.add(rbtnTTYes, "cell 2 8");
		
		rbtnTTNo = new JRadioButton("No");
		rbtnTTNo.setSelected(true);
		rbtnTTNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnTTNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				EnableTTValue(false);
			}
		});
		panel.add(rbtnTTNo, "cell 2 8");
		
		lblRatio = new JLabel("Ratio:");
		lblRatio.setEnabled(false);
		lblRatio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRatio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblRatio, "cell 2 8");
		
		txtTTRatio = new JTextField();
		txtTTRatio.setEnabled(false);
		txtTTRatio.setFont(new Font("Arial", Font.PLAIN, 13));
		txtTTRatio.setHorizontalAlignment(SwingConstants.CENTER);
		txtTTRatio.setText("2.0");
		txtTTRatio.setColumns(6);
		panel.add(txtTTRatio, "cell 2 8");
		
		lblSubRates = new JLabel("Set Substitution Rates:");
		lblSubRates.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSubRates.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSubRates, "flowx,cell 0 9 2 1,alignx right");
		
		btnSubsYes = new JRadioButton("Yes");
		btnSubsYes.setFont(new Font("Arial", Font.BOLD, 13));
		btnSubsYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CategoryToggle(true);
				CatsFileToggle(true);
				}
			});
		btnSubsYes.setSelected(false);
		panel.add(btnSubsYes, "cell 2 9");
		
		btnSubsNo = new JRadioButton("No");
		btnSubsNo.setFont(new Font("Arial", Font.BOLD, 13));
		btnSubsNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				CategoryToggle(false);
				CatsFileToggle(false);
				}
			});
		btnSubsNo.setSelected(true);
		panel.add(btnSubsNo, "cell 2 9");
		
		lblNumCats = new JLabel("Number of Substitution Categories:");
		lblNumCats.setEnabled(false);
		lblNumCats.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNumCats.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblNumCats, "flowx,cell 0 10 2 1,alignx right");

		cbxNumCats = new JComboBox();
		cbxNumCats.setEnabled(false);
		cbxNumCats.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cbxNumCats.setSelectedIndex(0);
		cbxNumCats.setFont(new Font("Arial", Font.PLAIN, 13));
		cbxNumCats.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				CategoryFileToggle();
			}
		});
		panel.add(cbxNumCats, "cell 2 10");
		
		lblCatNum = new JLabel("Category:");
		lblCatNum.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCatNum.setFont(new Font("Arial", Font.BOLD, 13));
		lblCatNum.setEnabled(false);
		panel.add(lblCatNum, "cell 2 10");
		
		cbxRateCatNum = new JComboBox();
		cbxRateCatNum.setModel(new DefaultComboBoxModel(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9"}));
		cbxRateCatNum.setFont(new Font("Arial", Font.PLAIN, 13));
		cbxRateCatNum.setEnabled(false);
		cbxRateCatNum.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				DisplayRateValue(Integer.parseInt(cbxRateCatNum.getSelectedItem().toString()));
			}
		});
		panel.add(cbxRateCatNum, "cell 2 10");
		
		lblRate = new JLabel("Rate:");
		lblRate.setHorizontalAlignment(SwingConstants.RIGHT);
		lblRate.setFont(new Font("Arial", Font.BOLD, 13));
		lblRate.setEnabled(false);
		panel.add(lblRate, "cell 2 10");
		
		txtRate = new JTextField();
		txtRate.setText("1.0");
		txtRate.setHorizontalAlignment(SwingConstants.CENTER);
		txtRate.setFont(new Font("Arial", Font.PLAIN, 13));
		txtRate.setEnabled(false);
		txtRate.setColumns(5);
		panel.add(txtRate, "cell 2 10");

		lblUseWeightsFor = new JLabel("Use weights for sites:");
		lblUseWeightsFor.setHorizontalAlignment(SwingConstants.RIGHT);
		lblUseWeightsFor.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseWeightsFor, "flowx,cell 0 11 2 1,alignx right");
		
		btnWeightYes = new JRadioButton("Yes");
		btnWeightYes.setFont(new Font("Arial", Font.BOLD, 13));
		btnWeightYes.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent arg0) {
			WeightToggle(true);
			}
		});
		btnWeightYes.setSelected(false);
		panel.add(btnWeightYes, "cell 2 11");
		
		btnWeightNo = new JRadioButton("No");
		btnWeightNo.setFont(new Font("Arial", Font.BOLD, 13));
		btnWeightNo.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent arg0) {
			WeightToggle(false);	
			}
		});
		btnWeightNo.setSelected(true);
		panel.add(btnWeightNo, "cell 2 11");
		
		lblUseEmpiricalBase = new JLabel("Use empirical base frequencies:");
		lblUseEmpiricalBase.setHorizontalAlignment(SwingConstants.RIGHT);
		lblUseEmpiricalBase.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblUseEmpiricalBase, "flowx,cell 0 12 2 1,alignx right");
		
		btnEmpYes = new JRadioButton("Yes");
		btnEmpYes.setFont(new Font("Arial", Font.BOLD, 13));
		btnEmpYes.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent arg0) {
			EmpToggle(true);				
			}
		});
		btnEmpYes.setSelected(true);
		panel.add(btnEmpYes, "cell 2 12");
		
		btnEmpNo = new JRadioButton("No");
		btnEmpNo.setFont(new Font("Arial", Font.BOLD, 13));
		btnEmpNo.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent arg0) {
			EmpToggle(false);				
			}
		});
		panel.add(btnEmpNo, "cell 2 12");
		
		lblBaseFreq = new JLabel("Base frequencies:");
		lblBaseFreq.setEnabled(false);
		lblBaseFreq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblBaseFreq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblBaseFreq, "flowx,cell 0 13 2 1,alignx right");
		
		lblBaseFreqA = new JLabel("A:");
		lblBaseFreqA.setEnabled(false);
		lblBaseFreqA.setFont(new Font("Arial", Font.BOLD, 13));
		lblBaseFreqA.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(lblBaseFreqA, "cell 2 13");
		
		txtBaseFreqA = new JTextField();
		txtBaseFreqA.setEnabled(false);
		txtBaseFreqA.setText("0.25");
		txtBaseFreqA.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqA.setBounds(292, 403, 43, 20);
		txtBaseFreqA.setColumns(6);
		panel.add(txtBaseFreqA, "cell 2 13");
	
		lblBaseFreqC = new JLabel("   C:  ");
		lblBaseFreqC.setEnabled(false);
		lblBaseFreqC.setFont(new Font("Arial", Font.BOLD, 13));
		lblBaseFreqC.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(lblBaseFreqC, "cell 2 13");
		
		txtBaseFreqC = new JTextField();
		txtBaseFreqC.setEnabled(false);
		txtBaseFreqC.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqC.setText("0.25");
		txtBaseFreqC.setColumns(6);
		panel.add(txtBaseFreqC, "cell 2 13");
		
		lblBaseFreqG = new JLabel("G:");
		lblBaseFreqG.setEnabled(false);
		lblBaseFreqG.setFont(new Font("Arial", Font.BOLD, 13));
		lblBaseFreqG.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(lblBaseFreqG, "cell 2 14");
		
		txtBaseFreqG = new JTextField();
		txtBaseFreqG.setEnabled(false);
		txtBaseFreqG.setText("0.25");
		txtBaseFreqG.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqG.setColumns(6);
		panel.add(txtBaseFreqG, "cell 2 14");
		
		lblBaseFreqTU = new JLabel("T/U:");
		lblBaseFreqTU.setEnabled(false);
		lblBaseFreqTU.setFont(new Font("Arial", Font.BOLD, 13));
		lblBaseFreqTU.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(lblBaseFreqTU, "cell 2 14");
		
		txtBaseFreqTU = new JTextField();
		txtBaseFreqTU.setEnabled(false);
		txtBaseFreqTU.setText("0.25");
		txtBaseFreqTU.setFont(new Font("Arial", Font.PLAIN, 13));
		txtBaseFreqTU.setColumns(6);
		panel.add(txtBaseFreqTU, "cell 2 14");

		lblDistMat = new JLabel("Distance matrix form:");
		lblDistMat.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDistMat.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDistMat, "flowx,cell 0 15 2 1,alignx right");
		
		cbxDistMat = new JComboBox();
		cbxDistMat.setFont(new Font("Arial", Font.PLAIN, 13));
		cbxDistMat.addItem("Square");
		cbxDistMat.addItem("Lower-triangle");
		cbxDistMat.addItem("Human-readable");
		panel.add(cbxDistMat, "cell 2 15");
		
		lblMultData = new JLabel("Analyze multiple data sets:");
		lblMultData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMultData, "flowx,cell 0 16 2 1,alignx right");
		
		btnMultYes = new JRadioButton("Yes");
		btnMultYes.setFont(new Font("Arial", Font.BOLD, 13));
		btnMultYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);				
				}
			});
		btnMultYes.setSelected(false);
		panel.add(btnMultYes, "cell 2 16");
		
		btnMultNo = new JRadioButton("No");
		btnMultNo.setFont(new Font("Arial", Font.BOLD, 13));
		btnMultNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);				
				}
			});
		btnMultNo.setSelected(true);
		panel.add(btnMultNo, "cell 2 16");
		
		lblKindData = new JLabel("Multiple data sets or multiple weights:");
		lblKindData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblKindData.setFont(new Font("Arial", Font.BOLD, 13));
		lblKindData.setEnabled(false);
		panel.add(lblKindData, "flowx,cell 0 17 2 1,alignx right");
		
		btnMultDataSet = new JRadioButton("Data sets");
		btnMultDataSet.setFont(new Font("Arial", Font.BOLD, 13));
		btnMultDataSet.setEnabled(false);
		btnMultDataSet.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(true);				
				}
			});
		btnMultDataSet.setSelected(true);
		panel.add(btnMultDataSet, "cell 2 17");
		
		btnMultWeight = new JRadioButton("Weights");
		btnMultWeight.setFont(new Font("Arial", Font.BOLD, 13));
		btnMultWeight.setEnabled(false);
		btnMultWeight.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DataWeightToggle(false);				
				}
			});
		panel.add(btnMultWeight, "cell 2 17");
		
		lblHowManyData = new JLabel("Number:");
		lblHowManyData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblHowManyData.setFont(new Font("Arial", Font.BOLD, 13));
		lblHowManyData.setEnabled(false);
		panel.add(lblHowManyData, "cell 2 18");
		
		txtDataWeight = new JTextField();
		txtDataWeight.setFont(new Font("Arial", Font.PLAIN, 13));
		txtDataWeight.setText("1");
		txtDataWeight.setEnabled(false);
		txtDataWeight.setColumns(5);
		panel.add(txtDataWeight, "cell 2 18");
		
		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setEnabled(false);
		lblInputSeq.setHorizontalAlignment(SwingConstants.RIGHT);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblInputSeq, "flowx,cell 0 19 2 1,alignx right");
		
		btnInputSeqYes = new JRadioButton("Interleaved");
		btnInputSeqYes.setEnabled(false);
		btnInputSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SeqToggle(true);				
				}
			});
		btnInputSeqYes.setSelected(true);
		panel.add(btnInputSeqYes, "cell 2 19");
		
		btnInputSeqNo = new JRadioButton("Sequential");
		btnInputSeqNo.setEnabled(false);
		btnInputSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		btnInputSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SeqToggle(false);				
				}
			});
		btnInputSeqNo.setSelected(false);
		panel.add(btnInputSeqNo, "cell 2 19");
		
		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintData, "flowx,cell 0 20 2 1,alignx right");
		
		btnPrintDataYes = new JRadioButton("Yes");
		btnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		btnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				doprintdataToggle(true);				
				}
			});
		btnPrintDataYes.setSelected(false);
		panel.add(btnPrintDataYes, "cell 2 20");
		
		btnPrintDataNo = new JRadioButton("No");
		btnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		btnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				doprintdataToggle(false);				
				}
			});
		btnPrintDataNo.setSelected(true);
		panel.add(btnPrintDataNo, "cell 2 20");
		
		lblDotDiff = new JLabel("Use dot-differencing to display them:");
		lblDotDiff.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDotDiff.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDotDiff, "flowx,cell 0 21 2 1,alignx right");
		
		btnDotDiffYes = new JRadioButton("Yes");
		btnDotDiffYes.setSelected(true);
		btnDotDiffYes.setHorizontalAlignment(SwingConstants.LEFT);
		btnDotDiffYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(true);
			}
		});
		btnDotDiffYes.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnDotDiffYes, "cell 2 21");
		
		btnDotDiffNo = new JRadioButton("No");
		btnDotDiffNo.setHorizontalAlignment(SwingConstants.LEFT);
		btnDotDiffNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DotDiffToggle(false);
			}
		});
		btnDotDiffNo.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnDotDiffNo, "cell 2 21");
		
		lblDoProgress = new JLabel("Display progress:");
		lblDoProgress.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDoProgress.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintData, "flowx,cell 0 22 2 1,alignx right");
		
		btnDoProgressYes = new JRadioButton("Yes");
		btnDoProgressYes.setFont(new Font("Arial", Font.BOLD, 13));
		btnDoProgressYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				doprintindToggle(true);				
				}
			});
		btnDoProgressYes.setSelected(true);
		panel.add(btnDoProgressYes, "cell 2 22");
		
		btnDoProgressNo = new JRadioButton("No");
		btnDoProgressNo.setFont(new Font("Arial", Font.BOLD, 13));
		btnDoProgressNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				doprintindToggle(false);				
				}
			});
		panel.add(btnDoProgressNo, "cell 2 22");
		
		btnExecute = new JButton("Execute");
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ratevalues[lastratecat] = Double.parseDouble(txtRate.getText());
				
				inputvals = new DnaDistData();
				
				inputvals.infile = txtInFile.getText();
				inputvals.outfile = txtOutFile.getText();
				inputvals.outfileopt = "w";
				inputvals.catfile = txtCatsFile.getText();
				inputvals.weightfile = txtWeightsFile.getText();
				inputvals.distanceoptions = cbxDistOpts.getSelectedItem().toString();
				inputvals.gammadistributed = GammaValue;
				inputvals.coefvar = Double.parseDouble(txtCoefVar.getText());
				inputvals.pctinvar = Double.parseDouble(txtPcntInvar.getText());
				inputvals.usettratio = rbtnTTYes.isSelected();
				inputvals.ttratio = Double.parseDouble(txtTTRatio.getText());
				inputvals.usesubstrates = btnSubsYes.isSelected();
				inputvals.numcats = Integer.parseInt(cbxNumCats.getSelectedItem().toString());
				inputvals.rate1 = ratevalues[0];
				inputvals.rate2 = ratevalues[1];
				inputvals.rate3 = ratevalues[2];
				inputvals.rate4 = ratevalues[3];
				inputvals.rate5 = ratevalues[4];
				inputvals.rate6 = ratevalues[5];
				inputvals.rate7 = ratevalues[6];
				inputvals.rate8 = ratevalues[7];
				inputvals.rate9 = ratevalues[8];
				inputvals.useweights = btnWeightYes.isSelected();
				inputvals.useEmpBF = btnEmpYes.isSelected();
				inputvals.basefreqA = Double.parseDouble(txtBaseFreqA.getText());
				inputvals.basefreqC = Double.parseDouble(txtBaseFreqC.getText());
				inputvals.basefreqG = Double.parseDouble(txtBaseFreqG.getText());
				inputvals.basefreqTU = Double.parseDouble(txtBaseFreqTU.getText());
				inputvals.distmatrix = cbxDistMat.getSelectedItem().toString();
				inputvals.dataset = Integer.parseInt(txtDataWeight.getText());
				inputvals.usemultdataset = btnMultYes.isSelected();
				inputvals.useData = btnMultDataSet.isSelected();
				inputvals.doseqint = btnInputSeqYes.isSelected();
				inputvals.doprintdata = btnPrintDataYes.isSelected();
				inputvals.dotdiff = btnDotDiffYes.isSelected();
				inputvals.doprintind = btnDoProgressYes.isSelected();
				
				btnExecute.setEnabled(false);	
				String title = "Dnadist Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread dnaDistThread = new Thread() {
						public void run() {
							runDnaDistThreads();
						}
			  	    };
			  	  dnaDistThread.start();
				}
				btnExecute.setEnabled(true);
			}
			
		});
		panel.add(btnExecute, "cell 2 23,alignx center");
		
		btnQuit = new JButton("Quit");
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmDnaDistControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		panel.add(btnQuit, "cell 2 23,alignx right");
	}
	
	public boolean checkInputVals(){
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
		
		if (inputvals.usesubstrates)
		{
			if (!test.FileAvailable(inputvals.catfile, "Categories"))
			{
				return false;
			}
		}
		
		if (inputvals.useweights)
		{
			if (!test.FileAvailable(inputvals.weightfile, "Weights"))
			{
				return false;
			}
		}

		if (inputvals.basefreqA < 0) {
			String msg1 = "Input value: Base Frequency for A cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.basefreqC < 0) {
			String msg1 = "Input value: Base Frequency for C cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.basefreqG < 0) {
			String msg1 = "Input value: Base Frequency for G cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		else if (inputvals.basefreqTU < 0) {
			String msg1 = "Input value: Base Frequency for T/U cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		double sum = inputvals.basefreqA + inputvals.basefreqC + inputvals.basefreqG + inputvals.basefreqTU;
		inputvals.basefreqA = inputvals.basefreqA / sum;
		inputvals.basefreqC = inputvals.basefreqC / sum;
		inputvals.basefreqG = inputvals.basefreqG / sum;
		inputvals.basefreqTU = inputvals.basefreqTU / sum;
		
		if (inputvals.coefvar <= 0){
			String msg1 = "Input value: Coefficient of Variation must be greater than 0.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if ((inputvals.pctinvar < 0) || (inputvals.pctinvar > 1)){
			String msg1 = "Input value: Fraction of Sites Invariant must be >= 0 and <= 1.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;		
		}
		
		if (inputvals.usesubstrates){
			// this could be compressed, but it keeps it in synch with the actual variables
			if ((inputvals.rate1 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 1 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
	
			if ((inputvals.numcats >= 2) && (inputvals.rate2 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 2 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
	
			if ((inputvals.numcats >= 3) && (inputvals.rate3 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 3 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
	
			if ((inputvals.numcats >= 4) && (inputvals.rate4 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 4 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
	
			if ((inputvals.numcats >= 5) && (inputvals.rate5 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 5 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
	
			if ((inputvals.numcats >= 6) && (inputvals.rate6 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 6 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
	
			if ((inputvals.numcats >= 7) && (inputvals.rate7 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 7 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
	
			if ((inputvals.numcats >= 8) && (inputvals.rate8 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 8 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
	
			if ((inputvals.numcats >= 9) && (inputvals.rate9 < 0)){
				String msg1 = "Input value: Substitution Rate for Category 9 must be >= 0.";
				JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
				return false;
			}
		}
		return true;
	}
	
	protected void runDnaDistThreads() {
    	try
    	{
    		// see if library exists
    		Native.loadLibrary("dnadist", DnaDist.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("DnaDist");
    		return;
    	}
		try 
		{
	  	    Thread dnaDistRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			DnaDist dnadist = (DnaDist) Native.loadLibrary("dnadist", DnaDist.class);
		  	        dnadist.dnadist(		
		  	        		inputvals.infile,
		  	        		inputvals.outfile,
		  	        		inputvals.outfileopt,
		  	        		inputvals.catfile,
		  	        		inputvals.weightfile,
		  	        		inputvals.distanceoptions,
		  	        		inputvals.gammadistributed,
		  	        		inputvals.coefvar,
		  	        		inputvals.pctinvar,
		  	        		inputvals.usettratio,
		  	        		inputvals.ttratio,
		  	        		inputvals.usesubstrates,
		  	        		inputvals.numcats,
		  	        		inputvals.rate1,
		  	        		inputvals.rate2,
		  	        		inputvals.rate3,
		  	        		inputvals.rate4,
		  	        		inputvals.rate5,
		  	        		inputvals.rate6,
		  	        		inputvals.rate7,
		  	        		inputvals.rate8,
		  	        		inputvals.rate9,
		  	        		inputvals.useweights,
		  	        		inputvals.useEmpBF,
		  	        		inputvals.basefreqA,
		  	        		inputvals.basefreqC,
		  	        		inputvals.basefreqG,
		  	        		inputvals.basefreqTU,
		  	        		inputvals.distmatrix,
		  	        		inputvals.usemultdataset,
		  	        		inputvals.useData,
		  	        		inputvals.dataset,
		  	        		inputvals.doseqint,
		  	        		inputvals.doprintdata,
		  	        		inputvals.dotdiff,
		  	        		inputvals.doprintind);		    		
		  	    };
	  	    };
	  	    dnaDistRunThread.start();

	  	    if (inputvals.doprintind)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (dnaDistRunThread.isAlive());
	  	    }
		} 
		catch (InterruptedException e) 
		{
			if (inputvals.doprintind)
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

