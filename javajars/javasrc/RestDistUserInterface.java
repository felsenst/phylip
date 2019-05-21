package phylip;
import java.awt.EventQueue;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import com.sun.jna.Library;
import com.sun.jna.Native;

import utilities.DisplayProgress;
import utilities.TestFileNames;

import java.awt.Font;
import java.awt.Color;

import net.miginfocom.swing.MigLayout;
import javax.swing.JScrollPane;
import java.awt.Dimension;
import javax.swing.JPanel;

public class RestDistUserInterface {
   public interface RestDist extends Library {
        public void restdist(
        		String infile,
        		String outfile,
        		String outfileopt,
        		boolean RestSites,
        		boolean OrigNeiLi,
        		boolean GammaDist,
        		double CoeffVar,
        		double TTRatio,
        		double SiteLen,
        		boolean SquareMat,
        		boolean AnalyzeMult,
        		int NumMult,
        		boolean Sequential,
        		boolean PrintData,
        		boolean PrintInd);
    }

	public class RestDistData {
		String infile;
		String outfile;
		String outfileopt;
		boolean RestSites;
		boolean OrigNeiLi;
		boolean GammaDist;
		double CoeffVar;
		double TTRatio;
		double SiteLen;
		boolean SquareMat;
		boolean AnalyzeMult;
		int NumMult;
		boolean Sequential;
		boolean PrintData;
		boolean PrintInd;
	}


	private RestDistData inputvals;
	private String inTitle;
	private String inCurdir;
	private DisplayProgress dp;
	private String filedir;
	private boolean phylipCall;

	private JFrame frmRestDistControls;
	private JLabel lblDistance;
	private JRadioButton rbtnRestSite;
	private JRadioButton rbtnFragment;
	private JLabel lblModel;
	private JRadioButton rbtnOriginal;
	private JRadioButton rbtnModified;
	private JLabel lblGamma;
	private JRadioButton rbtnGammaYes;
	private JRadioButton rbtnGammaNo;
	private JLabel lblInputSeq;
	private JRadioButton rbtnInputSeqNo;
	private JRadioButton rbtnInputSeqYes;
	private JLabel lblPrintData;
	private JRadioButton rbtnPrintDataYes;
	private JRadioButton rbtnPrintDataNo;
	private JLabel lblPrintInd;
	private JRadioButton rbtnPrintIndYes;
	private JRadioButton rbtnPrintIndNo;
	private JTextField txtInputFile;
	private JButton btnInputFile;
	private JButton btnExecute;
	private JButton btnQuit;
	private JTextField txtOutputFile;
	private JButton btnOutputFile;
	private JLabel lblTTRatio;
	private JTextField txtTTRatio;
	private JLabel lblMultData;
	private JRadioButton rbtnMultDataYes;
	private JRadioButton rbtnMultDataNo;
	private JLabel lblSiteLen;
	private JTextField txtSiteLen;
	private JLabel lblDMForm;
	private JRadioButton rbtnDMSquare;
	private JRadioButton rbtnDMLowerT;
	private JLabel lblNSets;
	private JTextField txtNSets;
	private JLabel lblCoeffVar;
	private JTextField txtCoeffVar;
	
	private JScrollPane scrollPane;
	private JPanel panel;


	/**
	 * Launch the application
	 */
	public static void main(final String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					RestDistUserInterface window = new RestDistUserInterface(args);
					window.frmRestDistControls.setVisible(true);
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

		int option = fileChooser.showOpenDialog(frmRestDistControls.getRootPane());
		if (option == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fileChooser.getSelectedFile();
			filedir = fileChooser.getCurrentDirectory().getAbsolutePath();
			file.setText(selectedFile.getPath());
		}
	}

	protected void SiteToggle(boolean isSite) {
		if (isSite) {
			rbtnRestSite.setSelected(true);
			rbtnFragment.setSelected(false);
		} else {
			rbtnRestSite.setSelected(false);
			rbtnFragment.setSelected(true);
		}
	}
	
	protected void ModelToggle(boolean isOriginal) {
		if (isOriginal) {
			rbtnOriginal.setSelected(true);
			rbtnModified.setSelected(false);
			lblGamma.setEnabled(false);
			rbtnGammaYes.setEnabled(false);
			rbtnGammaNo.setEnabled(false);
			lblCoeffVar.setEnabled(false);
			txtCoeffVar.setEnabled(false);
		} else {
			rbtnOriginal.setSelected(false);
			rbtnFragment.setSelected(true);
			lblGamma.setEnabled(true);
			rbtnGammaYes.setEnabled(true);
			rbtnGammaNo.setEnabled(true);
			lblCoeffVar.setEnabled(true);
			txtCoeffVar.setEnabled(true);
		}
	}

	protected void MultToggle(boolean isMult) {
		if (isMult) {
			rbtnMultDataYes.setSelected(true);
			rbtnMultDataNo.setSelected(false);
			lblInputSeq.setEnabled(true);
			rbtnInputSeqNo.setEnabled(true);
			rbtnInputSeqYes.setEnabled(true);
			lblNSets.setEnabled(true);
			txtNSets.setEnabled(true);
		} else {
			rbtnMultDataYes.setSelected(false);
			rbtnMultDataNo.setSelected(true);
			lblInputSeq.setEnabled(false);
			rbtnInputSeqNo.setEnabled(false);
			rbtnInputSeqYes.setEnabled(false);
			lblNSets.setEnabled(false);
			txtNSets.setEnabled(false);
		}
	}

	protected void GammaToggle(boolean isyes) {
		if (isyes) {
			rbtnGammaYes.setSelected(true);
			rbtnGammaNo.setSelected(false);
			lblCoeffVar.setEnabled(true);
			txtCoeffVar.setEnabled(true);
		} else {
			rbtnGammaYes.setSelected(false);
			rbtnGammaNo.setSelected(true);
			lblCoeffVar.setEnabled(false);
			txtCoeffVar.setEnabled(false);
		}
	}

	protected void DMSquare(boolean issquare) {
		if (issquare) {
			rbtnDMSquare.setSelected(true);
			rbtnDMLowerT.setSelected(false);
		} else {
			rbtnDMSquare.setSelected(false);
			rbtnDMLowerT.setSelected(true);
		}
	}

	protected void InputSeqToggle(boolean isInput) {
		if (isInput) {
			rbtnInputSeqNo.setSelected(true);
			rbtnInputSeqYes.setSelected(false);
		} else {
			rbtnInputSeqNo.setSelected(false);
			rbtnInputSeqYes.setSelected(true);
		}
	}

	protected void PrintDataToggle(boolean isPrintData) {
		if (isPrintData) {
			rbtnPrintDataYes.setSelected(true);
			rbtnPrintDataNo.setSelected(false);
		} else {
			rbtnPrintDataYes.setSelected(false);
			rbtnPrintDataNo.setSelected(true);
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
	
	/**
	 * Create the application.
	 */
	public RestDistUserInterface(String[] args) {
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
		
		filedir = System.getProperty("user.dir");

		frmRestDistControls = new JFrame();
		frmRestDistControls.setBackground(new Color(204, 255, 255));
		frmRestDistControls.setTitle("Restdist");
		frmRestDistControls.setBounds(100, 100, 630, 480);
		frmRestDistControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmRestDistControls.setPreferredSize(new Dimension(frmRestDistControls.getBounds().width, frmRestDistControls.getBounds().height));
		
		scrollPane = new JScrollPane();
		scrollPane.setPreferredSize(frmRestDistControls.getPreferredSize());
		frmRestDistControls.getContentPane().add(scrollPane);
		
		panel = new JPanel();
		panel.setPreferredSize(frmRestDistControls.getPreferredSize());
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
		panel.add(txtInputFile, "cell 1 0 2 1,growx");
	
		btnOutputFile = new JButton("Output File");
		btnOutputFile.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ChooseFile(txtOutputFile);
			}
		});
		btnOutputFile.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnOutputFile, "cell 0 1,growx");
		
		txtOutputFile = new JTextField();
		txtOutputFile.setFont(new Font("Arial", Font.PLAIN, 13));
		txtOutputFile.setText("outfile");
		panel.add(txtOutputFile, "cell 1 1 2 1,growx");

		lblDistance = new JLabel("Distances:");
		lblDistance.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDistance, "flowx,cell 0 2 2 1,alignx right");

		rbtnRestSite = new JRadioButton("Restriction sites");
		rbtnRestSite.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnRestSite.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnRestSite.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(true);
			}
		});
		rbtnRestSite.setSelected(true);
		panel.add(rbtnRestSite, "cell 2 2");

		rbtnFragment = new JRadioButton("Fragments");
		rbtnFragment.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnFragment.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnFragment.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				SiteToggle(false);
			}
		});
		panel.add(rbtnFragment, "cell 2 2");

		lblModel = new JLabel("Nei/Li model:");
		lblModel.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblModel, "flowx,cell 0 3 2 1,alignx right");

		rbtnOriginal = new JRadioButton("Original");
		rbtnOriginal.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnOriginal.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnOriginal.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ModelToggle(true);
			}
		});
		rbtnOriginal.setSelected(false);
		panel.add(rbtnOriginal, "cell 2 3");

		rbtnModified = new JRadioButton("Modified");
		rbtnModified.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnModified.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnModified.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				ModelToggle(false);
			}
		});
		rbtnModified.setSelected(true);
		panel.add(rbtnModified, "cell 2 3");

		lblGamma = new JLabel("Gamma distribution of rates among sites:");
		lblGamma.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblGamma, "flowx,cell 0 4 2 1,alignx right");

		rbtnGammaYes = new JRadioButton("Yes");
		rbtnGammaYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnGammaYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GammaToggle(true);
			}
		});
		panel.add(rbtnGammaYes, "cell 2 4");

		rbtnGammaNo = new JRadioButton("No");
		rbtnGammaNo.setSelected(true);
		rbtnGammaNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnGammaNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnGammaNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GammaToggle(false);
			}
		});
		panel.add(rbtnGammaNo, "cell 2 4");
		
		lblCoeffVar = new JLabel("Gamma coefficient of variation:");
		lblCoeffVar.setEnabled(false);
		lblCoeffVar.setHorizontalAlignment(SwingConstants.RIGHT);
		lblCoeffVar.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblCoeffVar, "flowx,cell 0 5 2 1,alignx right");

		txtCoeffVar = new JTextField();
		txtCoeffVar.setEnabled(false);
		txtCoeffVar.setText("1.0");
		txtCoeffVar.setFont(new Font("Arial", Font.PLAIN, 13));
		txtCoeffVar.setColumns(6);
		panel.add(txtCoeffVar, "cell 2 5");
		
		lblTTRatio = new JLabel("Transition/transversion ratio:");
		lblTTRatio.setHorizontalAlignment(SwingConstants.RIGHT);
		lblTTRatio.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblTTRatio, "flowx,cell 0 6 2 1,alignx right");
		
		txtTTRatio = new JTextField();
		txtTTRatio.setText("2.0");
		txtTTRatio.setFont(new Font("Arial", Font.PLAIN, 13));
		txtTTRatio.setColumns(6);
		panel.add(txtTTRatio, "cell 2 6");
		
		lblSiteLen = new JLabel("Site length:");
		lblSiteLen.setHorizontalAlignment(SwingConstants.RIGHT);
		lblSiteLen.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblSiteLen, "flowx,cell 0 7 2 1,alignx right");
		
		txtSiteLen = new JTextField();
		txtSiteLen.setText("6.0");
		txtSiteLen.setFont(new Font("Arial", Font.PLAIN, 13));
		txtSiteLen.setColumns(6);
		panel.add(txtSiteLen, "cell 2 7");
		
		lblDMForm = new JLabel("Distance Matrix Form:");
		lblDMForm.setHorizontalAlignment(SwingConstants.RIGHT);
		lblDMForm.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblDMForm, "flowx,cell 0 8 2 1,alignx right");
	
		rbtnDMSquare = new JRadioButton("Square");
		rbtnDMSquare.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DMSquare(true);
			}
		});
		rbtnDMSquare.setSelected(true);
		rbtnDMSquare.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDMSquare.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(rbtnDMSquare, "cell 2 8");
		
		rbtnDMLowerT = new JRadioButton("Lower Triangular");
		rbtnDMLowerT.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				DMSquare(false);
			}
		});
		rbtnDMLowerT.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnDMLowerT.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(rbtnDMLowerT, "cell 2 8");


		lblMultData = new JLabel("Analyze multiple data sets:");
		lblMultData.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblMultData, "flowx,cell 0 9 2 1,alignx right");

		rbtnMultDataYes = new JRadioButton("Yes");
		rbtnMultDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnMultDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(true);
			}
		});
		rbtnMultDataYes.setSelected(false);
		panel.add(rbtnMultDataYes, "cell 2 9");

		rbtnMultDataNo = new JRadioButton("No");
		rbtnMultDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnMultDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnMultDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				MultToggle(false);
			}
		});
		rbtnMultDataNo.setSelected(true);
		panel.add(rbtnMultDataNo, "cell 2 9");
		
		lblNSets = new JLabel("Number:\r\n");
		lblNSets.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNSets.setFont(new Font("Arial", Font.BOLD, 13));
		lblNSets.setEnabled(false);
		panel.add(lblNSets, "cell 2 9");
		
		txtNSets = new JTextField();
		txtNSets.setText("1");
		txtNSets.setFont(new Font("Arial", Font.PLAIN, 13));
		txtNSets.setEnabled(false);
		txtNSets.setColumns(6);
		panel.add(txtNSets, "cell 2 9");
	
		lblInputSeq = new JLabel("Input sequences:");
		lblInputSeq.setEnabled(false);
		lblInputSeq.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblInputSeq, "flowx,cell 0 10 2 1,alignx right");

		rbtnInputSeqNo = new JRadioButton("Interleaved");
		rbtnInputSeqNo.setEnabled(false);
		rbtnInputSeqNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnInputSeqNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnInputSeqNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(true);
			}
		});
		rbtnInputSeqNo.setSelected(true);
		panel.add(rbtnInputSeqNo, "cell 2 10");

		rbtnInputSeqYes = new JRadioButton("Sequential");
		rbtnInputSeqYes.setEnabled(false);
		rbtnInputSeqYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnInputSeqYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnInputSeqYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				InputSeqToggle(false);
			}
		});
		rbtnInputSeqYes.setSelected(false);
		panel.add(rbtnInputSeqYes, "cell 2 10");

		lblPrintData = new JLabel("Print out the data at start of run:");
		lblPrintData.setFont(new Font("Arial", Font.BOLD, 13));
		lblPrintData.setHorizontalAlignment(SwingConstants.RIGHT);
		panel.add(lblPrintData, "flowx,cell 0 11 2 1,alignx right");

		rbtnPrintDataYes = new JRadioButton("Yes");
		rbtnPrintDataYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintDataYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintDataYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(true);
			}
		});
		rbtnPrintDataYes.setSelected(false);
		panel.add(rbtnPrintDataYes, "cell 2 11");

		rbtnPrintDataNo = new JRadioButton("No");
		rbtnPrintDataNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintDataNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintDataNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintDataToggle(false);
			}
		});
		rbtnPrintDataNo.setSelected(true);
		panel.add(rbtnPrintDataNo, "cell 2 11");

		lblPrintInd = new JLabel("Display progress:");
		lblPrintInd.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(lblPrintInd, "flowx,cell 0 12 2 1,alignx right");

		rbtnPrintIndYes = new JRadioButton("Yes");
		rbtnPrintIndYes.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndYes.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintIndYes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(true);
			}
		});
		rbtnPrintIndYes.setSelected(true);
		panel.add(rbtnPrintIndYes, "cell 2 12");

		rbtnPrintIndNo = new JRadioButton("No");
		rbtnPrintIndNo.setFont(new Font("Arial", Font.BOLD, 13));
		rbtnPrintIndNo.setHorizontalAlignment(SwingConstants.LEFT);
		rbtnPrintIndNo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				PrintIndToggle(false);
			}
		});
		rbtnPrintIndNo.setSelected(false);
		panel.add(rbtnPrintIndNo, "cell 2 12");

		btnExecute = new JButton("Execute");
		btnExecute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				inputvals = new RestDistData();
				inputvals.infile = txtInputFile.getText();
				inputvals.outfile = txtOutputFile.getText();
				inputvals.outfileopt = "w";
				inputvals.RestSites = rbtnRestSite.isSelected();
				inputvals.OrigNeiLi = rbtnOriginal.isSelected();
				inputvals.GammaDist = rbtnGammaYes.isSelected();
				inputvals.CoeffVar =  Double.parseDouble(txtCoeffVar.getText());
				inputvals.TTRatio = Double.parseDouble(txtTTRatio.getText());
				inputvals.SiteLen = Double.parseDouble(txtSiteLen.getText());
				inputvals.SquareMat = rbtnDMSquare.isSelected();
				inputvals.AnalyzeMult = rbtnMultDataYes.isSelected();
				inputvals.NumMult =  Integer.parseInt(txtNSets.getText());;
				inputvals.Sequential = rbtnInputSeqYes.isSelected();
				inputvals.PrintData = rbtnPrintDataYes.isSelected();
				inputvals.PrintInd = rbtnPrintIndYes.isSelected();
				
				btnExecute.setEnabled(false);	
				String title = "Restdist Progress";
				String curdir = System.getProperty("user.dir");
				curdir += "/progress.txt";
				File fl = new File(curdir);
				fl.delete();
				inTitle = title;
				inCurdir = curdir;
	
				if (checkInputVals())
				{
					
			  	    Thread restdistThread = new Thread() {
						public void run() {
							runRestDistThreads();
						}
			  	    };
			  	  restdistThread.start();
				}
				btnExecute.setEnabled(true);
			}
		});
		btnExecute.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnExecute, "cell 2 13");

		btnQuit = new JButton("Quit");
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(phylipCall)
				{
					frmRestDistControls.dispose();
					return;
				}
				else
				{
					System.exit(0);
				}
			}
		});
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		panel.add(btnQuit, "cell 2 13");

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

		// check data
		if (inputvals.TTRatio < 0) {
			String msg1 = "Input value: Transition / Transversion Ratio cannot be negative.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		if (inputvals.CoeffVar <= 0.0) {
			String msg1 = "Input value: Coefficient of Variation must be greater than 0.";
			JOptionPane.showMessageDialog(null, msg1, "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}
	
	protected void runRestDistThreads() {
	   	try
	   	{
    		// see if library exists
    		Native.loadLibrary("restdist", RestDist.class);
		}
    	catch(UnsatisfiedLinkError e)
    	{
     		new TestFileNames().LibraryMissing("RestDist");
    		return;
    	}
		try 
		{
	  	    Thread restDistRunThread = new Thread() {
		  	      public void run() {
		  	    	  
		  			// at this point we hook into the C code
		  			RestDist restdist = (RestDist) Native.loadLibrary("restdist", RestDist.class);
		  			restdist.restdist(
		  	    		inputvals.infile,
		  	    		inputvals.outfile,
		  	    		inputvals.outfileopt,
		  	    		inputvals.RestSites,
		  	    		inputvals.OrigNeiLi,
		  	    		inputvals.GammaDist,
		  	    		inputvals.CoeffVar,
		  	    		inputvals.TTRatio,
		  	    		inputvals.SiteLen,
		  	    		inputvals.SquareMat,
		  	    		inputvals.AnalyzeMult,
		  	    		inputvals.NumMult,
		  	    		inputvals.Sequential,
		  	    		inputvals.PrintData,
		  	    		inputvals.PrintInd);
			  	    };
	  	    };
	  	    restDistRunThread.start();

	  	    if (inputvals.PrintInd)
	  	    {
		  	  	do
		  	  	{
					Thread.sleep(1000);
					updateProgress();
				} while (restDistRunThread.isAlive());
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
