package phylip;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
//import javax.swing.JScrollPane;

import org.eclipse.swt.*;
import org.eclipse.swt.widgets.*;
import org.eclipse.swt.layout.*;
//import net.miginfocom.swing.MigLayout;
import net.miginfocom.swing.MigLayout;

import clique.CliqueUserInterface;
import codml.CodMLUserInterface;
import consense.ConsenseUserInterface;
import contml.ContMLUserInterface;
import dnacomp.DnaCompUserInterface;
import dnadist.DnaDistUserInterface;
import dnainvar.DnaInvarUserInterface;
import dnaml.DnaMLUserInterface;
import dnamlk.DnaMLKUserInterface;
import dnapars.DnaParsUserInterface;
import dnapenny.DnaPennyUserInterface;
import dollop.DollopUserInterface;
import dolpenny.DolPennyUserInterface;
import drawgram.DrawgramUserInterface;
import drawtree.DrawtreeUserInterface;
import factor.FactorUserInterface;
import fitch.FitchUserInterface;
import gendist.GenDistUserInterface;
import kitsch.KitschUserInterface;
import mix.MixUserInterface;
import neighbor.NeighborUserInterface;
import pars.ParsUserInterface;
import penny.PennyUserInterface;
import proml.ProMLUserInterface;
import promlk.ProMLKUserInterface;
import protdist.ProtDistUserInterface;
import protpars.ProtParsUserInterface;
import restdist.RestDistUserInterface;
import retree.RetreeUserInterface;
import seqboot.SeqBootUserInterface;
import treedist.TreeDistUserInterface;

import com.sun.jna.Library;
import javax.swing.JButton;
import javax.swing.JTabbedPane;


public class PhylipUserInterface {
    public interface Phylip extends Library {
        public void  phylip(
		);
    }	
	private JFrame frmPhylipControls;
	private JTabbedPane tabbedPane;
	
	private JPanel pgmpanel;
	private JButton btnAboutPhylip;
	private JButton btnClique;
	private JButton btnCodml;
	private JButton btnConsense;
	private JButton btnContml;
	private JButton btnDnacomp;
	private JButton btnDnadist;
	private JButton btnDnainvar;
	private JButton btnDnaml;
	private JButton btnDnamlk;
	private JButton btnDnapars;
	private JButton btnDnapenny;
	private JButton btnDollop;
	private JButton btnDolpenny;
	private JButton btnDrawgram;
	private JButton btnDrawtree;
	private JButton btnFactor;
	private JButton btnFitch;
	private JButton btnGendist;
	private JButton btnKitsch;
	private JButton btnMix;
	private JButton btnNeighbor;
	private JButton btnPars;
	private JButton btnPenny;
	private JButton btnProml;
	private JButton btnPromlk;
	private JButton btnProtdist;
	private JButton btnProtpars;
	private JButton btnRestdist;
	private JButton btnRetree;
	private JButton btnSeqboot;
	private JButton btnTreedist;
	private JButton btnQuit;
	
	private JPanel pathpanel;
	private JButton btnBootstrap;
	private JLabel lblBootstrapPath;

	
	private JFrame frmAbout;
	private JButton btnAboutOK;


	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					PhylipUserInterface window = new PhylipUserInterface();
					window.frmPhylipControls.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public PhylipUserInterface() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		final String[] phylip = new String[] {"Phylip"};
		
		frmPhylipControls = new JFrame();
		frmPhylipControls.setBackground(new Color(204, 255, 255));
		frmPhylipControls.setFont(new Font("Arial", Font.BOLD, 13));
		frmPhylipControls.setTitle("Phylip");
		frmPhylipControls.setBounds(100, 100, 475, 450);
		frmPhylipControls.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmPhylipControls.setPreferredSize(new Dimension(frmPhylipControls.getBounds().width, frmPhylipControls.getBounds().height));
		
		tabbedPane = new JTabbedPane(JTabbedPane.TOP);
		tabbedPane.setPreferredSize(frmPhylipControls.getPreferredSize());

		frmPhylipControls.getContentPane().add(tabbedPane);
		
		pgmpanel = new JPanel();
		pgmpanel.setPreferredSize(new Dimension(455, 349));
		tabbedPane.addTab("Programs", pgmpanel);
		//scrollPane.setViewportView(panel);
		//scrollPane.setViewportView(panel);
		pgmpanel.setLayout(new MigLayout("", "[grow][][pref!,grow][pref!,grow]", "[grow][][][][][][][][][][]"));
		
		btnAboutPhylip = new JButton("About Phylip");
		btnAboutPhylip.setFont(new Font("Arial", Font.BOLD, 13));
		btnAboutPhylip.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
        		showAbout();
			}
		});
		
		pgmpanel.add(btnAboutPhylip, "cell 1 0 2 1,growx");
		
		btnClique = new JButton("Clique");
		btnClique.setFont(new Font("Arial", Font.PLAIN, 13));
		btnClique.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					CliqueOldUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnClique, "cell 0 1,growx");

		btnCodml = new JButton("CodML");
		btnCodml.setFont(new Font("Arial", Font.PLAIN, 13));
		btnCodml.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					CodMLUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnCodml, "cell 0 2,growx");
		
		btnConsense = new JButton("Consense");
		btnConsense.setFont(new Font("Arial", Font.PLAIN, 13));
		btnConsense.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					ConsenseUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnConsense, "cell 0 3,growx");
		
		btnContml = new JButton("Contml");
		btnContml.setFont(new Font("Arial", Font.PLAIN, 13));
		btnContml.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					ContMLUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnContml, "cell 0 4,growx");
		
		btnDnacomp = new JButton("Dnacomp");
		btnDnacomp.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDnacomp.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DnaCompUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDnacomp, "cell 0 5,growx");
		
		btnDnadist = new JButton("Dnadist");
		btnDnadist.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDnadist.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DnaDistUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDnadist, "cell 0 6,growx");
		
		btnDnainvar = new JButton("Dnainvar");
		btnDnainvar.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDnainvar.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DnaInvarUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDnainvar, "cell 0 7,growx");
		
		btnDnaml = new JButton("Dnaml");
		btnDnaml.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDnaml.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DnaMLUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDnaml, "cell 0 8,growx");
		
		btnDnamlk = new JButton("Dnamlk");
		btnDnamlk.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDnamlk.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DnaMLKUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDnamlk, "cell 1 1,growx");
		
		btnDnapars = new JButton("Dnapars");
		btnDnapars.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDnapars.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DnaParsUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDnapars, "cell 1 2,growx");
		
		btnDnapenny = new JButton("Dnapenny");
		btnDnapenny.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDnapenny.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DnaPennyUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDnapenny, "cell 1 3,growx");
		
		btnDollop = new JButton("Dollop");
		btnDollop.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDollop.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DollopUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDollop, "cell 1 4,growx");
		
		btnDolpenny = new JButton("Dolpenny");
		btnDolpenny.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDolpenny.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DolPennyUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDolpenny, "cell 1 5,growx");
		
		btnDrawgram = new JButton("Drawgram");
		btnDrawgram.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDrawgram.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DrawgramUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDrawgram, "cell 1 6,growx");
		
		btnDrawtree = new JButton("Drawtree");
		btnDrawtree.setFont(new Font("Arial", Font.PLAIN, 13));
		btnDrawtree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					DrawtreeUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnDrawtree, "cell 1 7,growx");
		
		btnFactor = new JButton("Factor");
		btnFactor.setFont(new Font("Arial", Font.PLAIN, 13));
		btnFactor.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					FactorUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnFactor, "cell 1 8,growx");
		
		btnFitch = new JButton("Fitch");
		btnFitch.setFont(new Font("Arial", Font.PLAIN, 13));
		btnFitch.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					FitchUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnFitch, "cell 2 1,growx");
		
		btnGendist = new JButton("Gendist");
		btnGendist.setFont(new Font("Arial", Font.PLAIN, 13));
		btnGendist.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					GenDistUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnGendist, "cell 2 2,growx");
		
		btnKitsch = new JButton("Kitsch");
		btnKitsch.setFont(new Font("Arial", Font.PLAIN, 13));
		btnKitsch.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					KitschUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnKitsch, "cell 2 3,growx");
		
		btnMix = new JButton("Mix");
		btnMix.setFont(new Font("Arial", Font.PLAIN, 13));
		btnMix.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					MixUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnMix, "cell 2 4,growx");
		
		btnNeighbor = new JButton("Neighbor");
		btnNeighbor.setFont(new Font("Arial", Font.PLAIN, 13));
		btnNeighbor.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					NeighborUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnNeighbor, "cell 2 5,growx");
		
		btnPars = new JButton("Pars");
		btnPars.setFont(new Font("Arial", Font.PLAIN, 13));
		btnPars.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					ParsUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnPars, "cell 2 6,growx");
		
		btnPenny = new JButton("Penny");
		btnPenny.setFont(new Font("Arial", Font.PLAIN, 13));
		btnPenny.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					PennyUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnPenny, "cell 2 7,growx");
		
		btnProml = new JButton("Proml");
		btnProml.setFont(new Font("Arial", Font.PLAIN, 13));
		btnProml.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					ProMLUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnProml, "cell 2 8,growx");
		
		btnPromlk = new JButton("Promlk");
		btnPromlk.setFont(new Font("Arial", Font.PLAIN, 13));
		btnPromlk.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					ProMLKUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnPromlk, "cell 3 1,growx");
		
		btnProtdist = new JButton("Protdist");
		btnProtdist.setFont(new Font("Arial", Font.PLAIN, 13));
		btnProtdist.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					ProtDistUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnProtdist, "cell 3 2,growx");
		
		btnProtpars = new JButton("Protpars");
		btnProtpars.setFont(new Font("Arial", Font.PLAIN, 13));
		btnProtpars.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					ProtParsUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnProtpars, "cell 3 3,growx");
		
		btnRestdist = new JButton("Restdist");
		btnRestdist.setFont(new Font("Arial", Font.PLAIN, 13));
		btnRestdist.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					RestDistUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnRestdist, "cell 3 4,growx");
		
		btnRetree = new JButton("Retree");
		btnRetree.setFont(new Font("Arial", Font.PLAIN, 13));
		btnRetree.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					RetreeUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnRetree, "cell 3 5,growx");
		
		btnSeqboot = new JButton("Seqboot");
		btnSeqboot.setFont(new Font("Arial", Font.PLAIN, 13));
		btnSeqboot.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					SeqBootUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnSeqboot, "cell 3 6,growx");
		
		btnTreedist = new JButton("Treedist");
		btnTreedist.setFont(new Font("Arial", Font.PLAIN, 13));
		btnTreedist.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					TreeDistUserInterface.main(phylip);
				} catch (Exception err) {
					err.printStackTrace();
				}
		 	}
		});
		pgmpanel.add(btnTreedist, "cell 3 7,growx");
		
		btnQuit = new JButton("Quit");
		btnQuit.setFont(new Font("Arial", Font.BOLD, 13));
		btnQuit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				System.exit(0);
		 	}
		});
		pgmpanel.add(btnQuit, "cell 3 9,growx");

	
		pathpanel = new JPanel();
		pathpanel.setPreferredSize(new Dimension(455, 349));
		tabbedPane.addTab("Pathways", pathpanel);
		pathpanel.setLayout(new MigLayout("", "[][grow]", "[][][][][][][][][][]"));
		
		
		btnBootstrap = new JButton("Bootstrap");
		btnBootstrap.setFont(new Font("Arial", Font.PLAIN, 13));
		btnBootstrap.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				final String[] bootstrap = new String[] {"Phylip", "bootstrap"};
				try {
					SeqBootUserInterface.main(bootstrap);
				} catch (Exception err) {
					err.printStackTrace();
				}
			}
		});
		pathpanel.add(btnBootstrap, "cell 0 0");

		lblBootstrapPath = new JLabel("Seqboot > Dnaml > Consense");
		lblBootstrapPath.setFont(new Font("Arial", Font.PLAIN, 13));
		pathpanel.add(lblBootstrapPath, "cell 1 0,growx");
	}
	
	public void showAbout() { 
		frmAbout = new JFrame();
		frmAbout.setVisible(true);
		frmAbout.setBounds(100, 50, 630, 300);
		frmAbout.setTitle("About Retree");
		frmAbout.getContentPane().setLayout(null);

		JLabel lblLine1 = new JLabel("Copyright 1993-2013. University of Washington and Joseph Felsenstein. All rights reserved.");
		lblLine1.setBounds(10, 10, 600, 14);
		frmAbout.getContentPane().add(lblLine1);
		
		JLabel lblLine2 = new JLabel("Permission is granted to reproduce, perform, and modify this program.");
		lblLine2.setBounds(10, 30, 600, 14);
		frmAbout.getContentPane().add(lblLine2);
   		
		JLabel lblLine3 = new JLabel("Permission is granted to distribute or provide access to this program provided that:");
		lblLine3.setBounds(10, 50, 600, 14);
		frmAbout.getContentPane().add(lblLine3);
   		
		JLabel lblLine4 = new JLabel("1) this copyright notice is not removed");
		lblLine4.setBounds(10,70, 600, 14);
		frmAbout.getContentPane().add(lblLine4);
   		
		JLabel lblLine5 = new JLabel("2) this program is not integrated with or called by any product or service that generates revenue");
		lblLine5.setBounds(10,90, 630, 14);
		frmAbout.getContentPane().add(lblLine5);
   		
		JLabel lblLine6 = new JLabel("3) your distribution of this program is free");
		lblLine6.setBounds(10,110, 600, 14);
		frmAbout.getContentPane().add(lblLine6);
   		
		JLabel lblLine7 = new JLabel("Any modified versions of this program that are distributed or accessible shall indicate");
		lblLine7.setBounds(10,150, 600, 14);
		frmAbout.getContentPane().add(lblLine7);
   		
		JLabel lblLine8 = new JLabel("that they are based on this program.  Educational institutions are granted permission");
		lblLine8.setBounds(10,170, 600, 14);
		frmAbout.getContentPane().add(lblLine8);
  		
		JLabel lblLine9 = new JLabel("to distribute this program to their students and staff for a fee to recover distribution costs.");
		lblLine9.setBounds(10,190, 600, 14);
		frmAbout.getContentPane().add(lblLine9);
  		
		JLabel lblLine10 = new JLabel("Permission requests for any other distribution of this program should be directed to:");
		lblLine10.setBounds(10,210, 600, 14);
		frmAbout.getContentPane().add(lblLine10);
  		
		JLabel lblLine11 = new JLabel("license (at) u.washington.edu.");
		lblLine11.setBounds(10,230, 600, 14);
		frmAbout.getContentPane().add(lblLine11);
		
		btnAboutOK = new JButton("OK");
		btnAboutOK.setBounds(500, 250, 84, 25);
		btnAboutOK.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				frmAbout.dispose();
				frmPhylipControls.repaint();        					
			}
		});
		frmAbout.getContentPane().add(btnAboutOK);
	}
}
