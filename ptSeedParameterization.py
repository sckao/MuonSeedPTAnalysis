import FWCore.ParameterSet.Config as cms

ptSeedParameterization = cms.PSet(
 
  vdouble CSC_01 = { 0.392612, 0.000000, 0.000000, 0.256925, 0.000000, 0.000000 } 
  vdouble CSC_12 = { -0.332655, 0.318077, 0.021135, -3.802048, 1.359701, 527.061915 } 
  vdouble CSC_02 = { -0.466204, 1.158276, -0.284877, 2.307452, 1.279234, -306.586626 } 
  vdouble CSC_13 = { 0.161190, -0.433689, 0.321022, -2.464092, 1.622968, 360.885026 } 
  vdouble CSC_03 = { -0.256329, 0.938594, -0.201104, 2.446152, 1.457664, -323.756776 } 
  vdouble CSC_14 = { -0.600598, 1.243747, -0.262080, 3.356669, 1.821404, -452.853439 } 
  vdouble CSC_23 = { -0.088570, 0.067946, 0.013028, 214.374032, -5115.096735, 30640.002899 } 
  vdouble CSC_24 = { 0.486351, -0.506147, 0.161176, -1.673370, 5.273864, 377.480225 } 
  vdouble CSC_34 = { 0.678265, -0.676532, 0.172211, -3491.527274, 84116.760862, -505816.660089 } 
  vdouble DT_12 = { 0.182946, 0.054820, -0.085796, 12.546145, -256.056605, 1321.642971 } 
  vdouble DT_13 = { 0.315003, 0.067256, -0.125512, 0.186846, 0.539971, -8.785914 } 
  vdouble DT_14 = { 0.361867, 0.035392, -0.091604, 0.292670, 0.688432, -17.170856 } 
  vdouble DT_23 = { 0.129178, 0.026056, -0.059270, 27.602430, -559.322097, 2854.990753 } 
  vdouble DT_24 = { 0.173530, 0.039266, -0.088289, 0.660922, 1.104137, -48.399435 } 
  vdouble DT_34 = { 0.039513, 0.017686, -0.025632, 51.966228, -1058.615873, 5466.806954 } 
  vdouble OL_1213 = { 0.206344, 0.000000, 0.000000, 0.210755, 0.000000, 0.000000 } 
  vdouble OL_1222 = { 0.993100, -0.664834, 0.000000, -1.322209, 17.769488, 0.000000 } 
  vdouble OL_1232 = { 0.183154, 0.000000, 0.000000, 0.511705, -1.689062, 0.000000 } 
  vdouble OL_2213 = { 0.100036, 0.000000, 0.000000, 0.324587, 0.000000, 0.000000 } 
  vdouble OL_2222 = { 0.107933, 0.000000, 0.000000, 0.573849, 0.000000, 0.000000 } 
  vdouble SME_11 = { 3.089745, -1.206591, 0.032380, 244.655663, -5718.138491, 33462.543590 } 
  vdouble SME_12 = { 0.060489, 0.585768, 0.000000, 0.405479, 0.150504, 0.000000 } 
  vdouble SME_13 = { -1.453585, 1.713315, 0.000000, 6.234882, -59.697311, 0.000000 } 
  vdouble SME_21 = { -0.374400, 0.982096, -0.284961, 332.538723, -7822.926723, 46219.959517 } 
  vdouble SME_22 = { -1.894764, 2.391133, -0.585444, 6127.897542, -140666.238615, 807384.478412 } 
  vdouble SME_31 = { -0.674846, 0.556384, -0.086880, -937.952107, 21155.546349, -117390.639201 } 
  vdouble SME_32 = { -0.441514, 0.605367, -0.187861, 19849.178146, -459240.087996, 2658160.556684 } 
  vdouble SME_41 = { -2.618203, 2.373979, -0.525277, 35123.036510, -861695.545076, 5287778.986458 } 
  vdouble SMB_10 = { 1.383941, -0.025433, 0.000000, -0.109540, 2.498832, 0.000000 } 
  vdouble SMB_11 = { 1.237265, 0.756238, -0.833452, 0.074032, 0.452398, 1.494974 } 
  vdouble SMB_12 = { 2.188682, -1.017835, 0.000000, -3.264964, 37.837152, 0.000000 } 
  vdouble SMB_20 = { 1.012883, -0.070366, 0.000000, 0.177533, 0.057712, 0.000000 } 
  vdouble SMB_21 = { 1.053309, -0.140623, 0.000000, 0.463200, -2.967335, 0.000000 } 
  vdouble SMB_22 = { 1.454569, -0.734477, 0.000000, 2.346655, -23.065756, 0.000000 } 
  vdouble SMB_30 = { 0.507874, -0.023191, 0.000000, 2.238747, -18.314311, 0.000000 } 
  vdouble SMB_31 = { 0.553217, -0.159012, 0.000000, 2.802067, -24.749146, 0.000000 } 
  vdouble SMB_32 = { 0.671291, -0.328846, 0.000000, -1.503412, 21.832958, 0.000000 } 
)