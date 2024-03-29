import FWCore.ParameterSet.Config as cms

dphiScale = cms.PSet(
    # ============================================================================== 
    #    idx                         p0        p1        p2        r1         r2 
    # ------------------------------------------------------------------------------ 

  CSC_01_1_scale = cms.vdouble( 0.386926, -1.180591, 0.000000, -3.051205, 0.000000 ), 
  CSC_12_1_scale = cms.vdouble( 0.144868, -4.216645, 0.000000, -29.106790, 0.000000 ), 
  CSC_12_2_scale = cms.vdouble( 0.716967, -1.340694, 0.000000, -1.869953, 0.000000 ), 
  CSC_12_3_scale = cms.vdouble( 0.716967, -1.340694, 0.000000, -1.869953, 0.000000 ), 
  CSC_13_2_scale = cms.vdouble( 0.225079, -4.871236, 0.000000, -21.642340, 0.000000 ), 
  CSC_13_3_scale = cms.vdouble( 0.840912, -1.535245, 0.000000, -1.825690, 0.000000 ), 
  CSC_14_3_scale = cms.vdouble( 0.872768, -2.054468, 0.000000, -2.353969, 0.000000 ), 
  CSC_23_1_scale = cms.vdouble( 0.015700, -0.000100, 0.000000, -0.006371, 0.000000 ), 
  CSC_23_2_scale = cms.vdouble( 0.105480, -3.296721, 0.000000, -31.254558, 0.000000 ), 
  CSC_24_1_scale = cms.vdouble( 0.137168, -4.961803, 0.000000, -36.173204, 0.000000 ), 
  CSC_34_1_scale = cms.vdouble( 0.013730, -10.000000, 0.000000, -728.353036, 0.000000 ), 
  OL_1213_0_scale = cms.vdouble( 0.292489, -3.652815, 0.000000, -12.488734, 0.000000 ), 
  OL_1222_0_scale = cms.vdouble( 0.241839, -4.893208, 0.000000, -20.233327, 0.000000 ), 
  OL_1232_0_scale = cms.vdouble( 0.198508, -4.860602, 0.000000, -24.485651, 0.000000 ), 
  OL_2213_0_scale = cms.vdouble( 0.124530, -6.169639, 0.000000, -49.543234, 0.000000 ), 
  OL_2222_0_scale = cms.vdouble( 0.119692, -6.783315, 0.000000, -56.673236, 0.000000 ), 
  DT_12_1_scale = cms.vdouble( 0.197696, -2.591764, 0.000000, -13.109831, 0.000000 ), 
  DT_12_2_scale = cms.vdouble( 0.173231, -2.427548, 0.000000, -14.013322, 0.000000 ), 
  DT_13_1_scale = cms.vdouble( 0.340967, -3.727984, 0.000000, -10.933548, 0.000000 ), 
  DT_13_2_scale = cms.vdouble( 0.322764, -3.440971, 0.000000, -10.660939, 0.000000 ), 
  DT_14_1_scale = cms.vdouble( 0.393400, -5.128164, 0.000000, -13.035502, 0.000000 ), 
  DT_14_2_scale = cms.vdouble( 0.378315, -4.165588, 0.000000, -11.010894, 0.000000 ), 
  DT_23_1_scale = cms.vdouble( 0.141913, -4.008556, 0.000000, -28.246634, 0.000000 ), 
  DT_23_2_scale = cms.vdouble( 0.121109, -3.724435, 0.000000, -30.752771, 0.000000 ), 
  DT_24_1_scale = cms.vdouble( 0.196086, -6.828230, 0.000000, -34.822666, 0.000000 ), 
  DT_24_2_scale = cms.vdouble( 0.184759, -5.646392, 0.000000, -30.560800, 0.000000 ), 
  DT_34_1_scale = cms.vdouble( 0.050022, -10.000000, 0.000000, -199.910415, 0.000000 ), 
  DT_34_2_scale = cms.vdouble( 0.048093, -10.000000, 0.000000, -207.928367, 0.000000 ), 
  SMB_10_0_scale = cms.vdouble( 1.311124, 2.684074, 0.000000, 2.047154, 0.000000 ), 
  SMB_11_0_scale = cms.vdouble( 1.332057, 2.834192, 0.000000, 2.127682, 0.000000 ), 
  SMB_12_0_scale = cms.vdouble( 1.202334, 2.381430, 0.000000, 1.980673, 0.000000 ), 
  SMB_20_0_scale = cms.vdouble( 0.976241, 1.412898, 0.000000, 1.447285, 0.000000 ), 
  SMB_21_0_scale = cms.vdouble( 0.957840, 1.487004, 0.000000, 1.552456, 0.000000 ), 
  SMB_22_0_scale = cms.vdouble( 0.860039, 1.208654, 0.000000, 1.405348, 0.000000 ), 
  SMB_30_0_scale = cms.vdouble( 0.542504, -3.275737, 0.000000, -6.038180, 0.000000 ), 
  SMB_31_0_scale = cms.vdouble( 0.524577, -2.915318, 0.000000, -5.557462, 0.000000 ), 
  SMB_32_0_scale = cms.vdouble( 0.466847, -2.581907, 0.000000, -5.530521, 0.000000 ), 
  SME_11_0_scale = cms.vdouble( 0.867166, 1.184825, 0.000000, 1.366319, 0.000000 ), 
  SME_12_0_scale = cms.vdouble( 0.873350, 1.967434, 0.000000, 2.252745, 0.000000 ), 
  SME_13_0_scale = cms.vdouble( 0.553659, 0.057663, 0.000000, 0.104148, 0.000000 ), 
  SME_21_0_scale = cms.vdouble( 0.447316, -0.100830, 0.000000, -0.225411, 0.000000 ), 
  SME_22_0_scale = cms.vdouble( 0.308319, -2.469849, 0.000000, -8.010697, 0.000000 )
)
