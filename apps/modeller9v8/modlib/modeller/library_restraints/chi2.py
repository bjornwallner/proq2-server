#   residue   spline values
chi2 = [
    ( 'ASP', (1.8744, 1.9556, 2.0399, 2.1535, 2.3003, 2.4334, 2.5631, 2.6677, 2.7641, 2.8517, 2.9476, 3.0379, 3.1161, 3.2110, 3.2797, 3.3663, 3.4470, 3.4633, 3.4798, 3.4556, 3.4193, 3.4183, 3.3853, 3.3524, 3.3225, 3.2571, 3.2050, 3.1499, 3.0803, 3.0437, 2.9979, 2.9721, 2.9503, 2.8552, 2.7288, 2.5424, 2.3376, 2.1762, 2.0394, 1.9338, 1.8382, 1.7242, 1.6007, 1.4657, 1.3338, 1.2124, 1.0968, 0.9910, 0.8967, 0.8161, 0.7475, 0.6899, 0.6404, 0.5910, 0.5453, 0.5014, 0.4552, 0.4070, 0.3527, 0.2913, 0.2245, 0.1593, 0.1008, 0.0525, 0.0184, 0.0000, 0.0022, 0.0297, 0.0855, 0.1653, 0.2628, 0.3687, 0.4737, 0.5805, 0.6774, 0.7656, 0.8440, 0.9064, 0.9712, 1.0336, 1.1021, 1.1741, 1.2432, 1.3149, 1.3823, 1.4524, 1.5191, 1.5825, 1.6267, 1.6583, 1.6690, 1.6497, 1.6297, 1.5953, 1.5656, 1.5584, 1.5668, 1.5946, 1.6482, 1.7192, 1.7978, 1.8814, 1.9579, 2.0276, 2.0980, 2.2022, 2.3405, 2.5077, 2.7132, 2.8794, 2.9924, 3.0626, 3.0893, 3.1145, 3.1356, 3.1515, 3.1579, 3.1579, 3.1571, 3.1356, 3.0954, 3.0670, 3.0223, 2.9993, 2.9775, 2.9340, 2.8911, 2.8276, 2.7554, 2.6603, 2.5631, 2.4520, 2.3288, 2.2072, 2.0742, 1.9372, 1.8122, 1.7028, 1.6085, 1.5471, 1.5103, 1.4901, 1.4931, 1.4991) ),
    ( 'GLU', (0.0821, 0.1640, 0.2776, 0.4316, 0.6551, 0.9080, 1.1844, 1.4829, 1.7815, 2.0874, 2.3819, 2.6447, 2.8946, 3.1001, 3.2941, 3.4809, 3.6287, 3.7758, 3.8919, 3.9541, 4.0133, 4.0725, 4.1162, 4.1641, 4.1915, 4.1558, 4.0917, 4.0017, 3.8738, 3.7174, 3.5305, 3.3126, 3.0851, 2.8634, 2.6687, 2.5022, 2.3590, 2.2383, 2.1167, 1.9955, 1.8714, 1.7439, 1.6291, 1.5250, 1.4483, 1.4004, 1.3837, 1.4098, 1.4777, 1.6032, 1.7922, 2.0497, 2.3732, 2.7434, 3.1439, 3.5350, 3.9250, 4.2753, 4.5858, 4.8458, 5.0261, 5.2099, 5.3464, 5.4465, 5.5125, 5.5493, 5.5789, 5.6228, 5.6640, 5.6640, 5.6547, 5.6409, 5.5962, 5.6006, 5.6138, 5.6318, 5.6781, 5.6593, 5.6363, 5.5875, 5.4966, 5.4578, 5.3635, 5.2369, 5.1219, 4.9274, 4.7074, 4.4605, 4.1242, 3.7506, 3.3675, 2.9764, 2.6199, 2.3121, 2.0608, 1.8774, 1.7590, 1.6986, 1.6874, 1.7153, 1.7702, 1.8443, 1.9187, 1.9877, 2.0494, 2.0993, 2.1720, 2.2671, 2.3932, 2.5714, 2.7909, 3.0482, 3.3383, 3.6238, 3.8408, 4.0026, 4.1093, 4.1830, 4.2466, 4.2871, 4.3073, 4.2954, 4.2638, 4.2132, 4.1343, 4.0017, 3.8761, 3.7107, 3.5081, 3.2976, 3.0417, 2.7744, 2.4886, 2.1835, 1.8682, 1.5478, 1.2400, 0.9528, 0.6905, 0.4629, 0.2752, 0.1459, 0.0574, 0.0000) ),
    ( 'PHE', (3.9480, 4.0061, 4.0587, 4.1107, 4.1620, 4.2183, 4.2659, 4.3005, 4.3132, 4.3081, 4.3030, 4.2880, 4.2930, 4.3106, 4.3183, 4.3286, 4.3209, 4.2683, 4.1774, 4.0779, 3.9547, 3.8246, 3.7136, 3.5837, 3.4482, 3.3156, 3.1672, 3.0024, 2.8326, 2.6492, 2.4406, 2.2347, 1.9977, 1.6829, 1.3613, 1.0558, 0.8042, 0.6438, 0.5652, 0.5573, 0.6090, 0.7028, 0.8068, 0.9291, 1.0593, 1.1867, 1.3394, 1.4892, 1.6443, 1.8004, 1.9514, 2.0997, 2.2344, 2.3618, 2.4556, 2.5325, 2.5955, 2.6381, 2.6716, 2.6876, 2.6931, 2.6966, 2.7058, 2.7213, 2.7615, 2.8091, 2.8549, 2.9275, 2.9820, 3.0388, 3.1044, 3.1345, 3.1584, 3.1836, 3.2070, 3.2589, 3.3109, 3.3540, 3.3777, 3.3423, 3.2951, 3.2291, 3.1560, 3.1013, 3.0395, 2.9554, 2.8603, 2.7508, 2.6164, 2.4806, 2.3259, 2.1526, 1.9782, 1.7883, 1.5908, 1.3853, 1.1703, 0.9612, 0.7588, 0.5726, 0.4088, 0.2637, 0.1468, 0.0561, 0.0015, 0.0000, 0.0522, 0.1697, 0.3469, 0.5554, 0.7772, 0.9884, 1.1711, 1.3464, 1.5221, 1.6809, 1.8416, 1.9889, 2.1259, 2.2708, 2.4131, 2.5596, 2.7032, 2.8367, 2.9314, 3.0051, 3.0722, 3.1221, 3.1977, 3.2678, 3.3044, 3.3289, 3.3433, 3.3569, 3.3747, 3.3979, 3.4020, 3.4020, 3.4301, 3.4787, 3.5463, 3.6310, 3.7052, 3.7647) ),
    ( 'HIS', (1.5542, 1.6170, 1.6661, 1.6960, 1.7546, 1.7972, 1.8654, 1.9654, 2.0429, 2.1235, 2.1893, 2.2374, 2.2785, 2.3114, 2.3475, 2.3681, 2.3828, 2.3681, 2.2833, 2.1711, 2.0481, 1.9339, 1.8461, 1.7642, 1.6776, 1.5744, 1.4489, 1.3131, 1.1754, 1.0336, 0.9139, 0.8014, 0.6935, 0.5997, 0.4944, 0.3936, 0.3016, 0.2156, 0.1425, 0.0834, 0.0367, 0.0065, 0.0000, 0.0146, 0.0582, 0.1269, 0.2266, 0.3556, 0.5087, 0.6955, 0.8986, 1.1207, 1.3523, 1.5712, 1.7677, 1.9475, 2.1033, 2.2550, 2.4042, 2.5136, 2.5939, 2.6516, 2.6926, 2.7337, 2.7893, 2.8216, 2.8532, 2.8876, 2.9071, 2.9379, 2.9793, 3.0164, 3.0866, 3.1620, 3.2189, 3.2612, 3.2587, 3.2336, 3.1972, 3.1574, 3.1258, 3.0951, 3.0612, 3.0406, 3.0046, 2.9324, 2.8265, 2.6489, 2.4315, 2.2229, 2.0037, 1.7884, 1.5773, 1.3599, 1.1489, 0.9493, 0.7721, 0.6120, 0.4845, 0.3913, 0.3262, 0.2939, 0.2806, 0.2844, 0.3040, 0.3460, 0.4133, 0.5061, 0.6147, 0.7212, 0.8188, 0.8991, 0.9764, 1.0581, 1.1492, 1.2455, 1.3535, 1.4665, 1.5815, 1.7040, 1.8008, 1.8743, 1.9265, 1.9640, 2.0045, 2.0459, 2.0796, 2.1009, 2.1235, 2.1400, 2.1450, 2.1178, 2.0542, 1.9696, 1.8724, 1.7949, 1.7148, 1.6412, 1.5806, 1.5302, 1.5075, 1.5010, 1.5110, 1.5208) ),
    ( 'ILE', (1.3757, 1.6187, 1.9148, 2.2781, 2.8415, 3.3676, 3.8277, 4.2561, 4.5721, 4.8221, 5.0022, 5.1193, 5.2369, 5.3316, 5.4283, 5.5267, 5.6264, 5.6852, 5.7335, 5.7695, 5.7768, 5.8144, 5.8260, 5.8182, 5.7992, 5.7768, 5.7658, 5.7695, 5.7549, 5.6920, 5.6044, 5.4520, 5.3082, 5.1744, 5.0492, 4.9380, 4.8066, 4.6150, 4.3463, 4.0220, 3.6431, 3.2462, 2.8474, 2.4555, 2.0968, 1.8031, 1.5963, 1.4930, 1.5065, 1.6381, 1.8875, 2.2462, 2.6900, 3.1855, 3.6750, 4.1283, 4.5235, 4.8321, 5.0818, 5.2764, 5.4257, 5.5470, 5.6586, 5.7229, 5.7731, 5.8377, 5.8736, 5.9236, 5.9944, 6.0318, 6.0657, 6.0806, 6.0608, 6.0559, 6.0608, 6.0559, 6.0608, 6.0657, 6.0462, 6.0271, 5.9990, 5.9540, 5.9452, 5.9452, 5.9109, 5.8615, 5.7658, 5.6359, 5.4928, 5.2922, 4.9988, 4.6551, 4.2666, 3.8718, 3.5359, 3.2613, 3.0737, 2.9846, 2.9740, 3.0397, 3.1653, 3.3152, 3.4735, 3.6253, 3.7422, 3.8351, 3.9167, 3.9716, 4.0303, 4.1072, 4.1883, 4.2981, 4.4148, 4.5099, 4.6024, 4.6782, 4.7602, 4.8539, 4.9316, 4.9787, 4.9476, 4.8408, 4.7184, 4.5888, 4.4592, 4.3384, 4.1913, 4.0062, 3.7835, 3.5184, 3.1889, 2.8048, 2.3749, 1.9292, 1.4936, 1.0904, 0.7389, 0.4471, 0.2229, 0.0752, 0.0113, 0.0000, 0.0357, 0.1057) ),
    ( 'LYS', (0.0585, 0.1561, 0.2932, 0.4755, 0.7347, 1.0167, 1.3067, 1.5918, 1.8547, 2.1132, 2.3547, 2.5796, 2.7915, 2.9621, 3.1053, 3.2307, 3.3520, 3.4773, 3.5844, 3.6844, 3.7658, 3.8524, 3.9708, 4.0598, 4.1189, 4.1208, 4.0756, 4.0417, 4.0123, 3.9788, 3.9166, 3.8143, 3.6952, 3.5704, 3.4519, 3.3435, 3.2123, 3.0782, 2.9181, 2.7237, 2.5322, 2.3134, 2.1041, 1.9309, 1.7810, 1.6961, 1.6743, 1.7077, 1.8247, 2.0036, 2.2441, 2.5569, 2.8994, 3.2615, 3.6151, 3.9334, 4.2016, 4.4429, 4.6680, 4.8458, 5.0205, 5.1751, 5.2992, 5.3969, 5.4761, 5.5314, 5.5621, 5.5779, 5.6143, 5.6477, 5.6736, 5.7092, 5.6912, 5.6824, 5.6649, 5.6225, 5.6184, 5.6061, 5.6061, 5.6225, 5.5979, 5.5466, 5.4869, 5.4102, 5.2962, 5.1645, 5.0070, 4.8195, 4.5971, 4.3263, 4.0298, 3.6850, 3.3295, 3.0032, 2.7064, 2.4874, 2.3527, 2.2884, 2.3084, 2.3857, 2.5041, 2.6642, 2.8314, 3.0185, 3.2123, 3.3882, 3.5678, 3.7055, 3.8042, 3.8928, 3.9598, 4.0206, 4.0809, 4.1171, 4.1282, 4.1328, 4.1470, 4.1681, 4.1906, 4.2259, 4.2106, 4.1661, 4.1034, 3.9999, 3.8921, 3.7665, 3.6263, 3.4730, 3.3076, 3.1407, 2.9592, 2.7638, 2.5633, 2.3256, 2.0688, 1.7949, 1.4987, 1.2060, 0.9134, 0.6368, 0.3926, 0.2182, 0.0909, 0.0000) ),
    ( 'LEU', (0.5398, 0.7126, 0.9355, 1.2271, 1.6572, 2.1055, 2.5440, 2.9208, 3.2067, 3.4091, 3.5542, 3.6866, 3.8104, 3.9482, 4.0776, 4.2268, 4.3910, 4.5784, 4.8104, 5.0035, 5.2000, 5.3623, 5.4754, 5.5979, 5.6600, 5.7092, 5.7609, 5.7878, 5.8186, 5.8123, 5.8000, 5.7698, 5.7120, 5.6156, 5.4622, 5.3084, 5.1428, 4.9911, 4.8265, 4.6303, 4.4299, 4.2332, 4.0792, 3.9765, 3.9315, 3.9453, 3.9925, 4.0721, 4.1687, 4.2585, 4.3614, 4.4616, 4.5569, 4.6722, 4.8036, 4.9549, 5.1239, 5.2990, 5.4257, 5.5208, 5.5781, 5.5880, 5.5954, 5.5855, 5.5514, 5.5208, 5.4843, 5.4194, 5.3564, 5.2715, 5.1539, 5.0372, 4.8788, 4.7168, 4.5507, 4.3850, 4.2389, 4.0958, 3.9785, 3.8700, 3.7786, 3.7008, 3.6276, 3.5601, 3.5008, 3.4463, 3.3708, 3.2674, 3.0846, 2.8023, 2.4229, 1.9946, 1.5765, 1.2166, 0.9415, 0.7612, 0.6825, 0.7017, 0.8169, 1.0196, 1.2973, 1.6383, 2.0134, 2.4013, 2.7693, 3.0820, 3.3521, 3.5886, 3.8125, 4.0473, 4.2717, 4.4786, 4.6652, 4.8127, 4.9457, 5.0430, 5.1023, 5.1636, 5.2034, 5.2256, 5.2430, 5.2256, 5.1604, 5.0444, 4.8643, 4.6702, 4.4362, 4.1992, 3.9550, 3.6792, 3.4119, 3.1197, 2.7821, 2.4045, 1.9853, 1.5565, 1.1580, 0.8066, 0.5161, 0.2956, 0.1446, 0.0567, 0.0131, 0.0000) ),
    ( 'MET', (0.0274, 0.0911, 0.1911, 0.3383, 0.5603, 0.8402, 1.1650, 1.5170, 1.8799, 2.2131, 2.5069, 2.7653, 2.9795, 3.1860, 3.3722, 3.5493, 3.7181, 3.8387, 3.9341, 3.9921, 4.0030, 4.0225, 4.0567, 4.0861, 4.1041, 4.0891, 4.0337, 3.9680, 3.9038, 3.8481, 3.7999, 3.7038, 3.5687, 3.3693, 3.1248, 2.8729, 2.6418, 2.4443, 2.2741, 2.1262, 1.9532, 1.7555, 1.5343, 1.3012, 1.0975, 0.9263, 0.8064, 0.7534, 0.7626, 0.8504, 1.0173, 1.2511, 1.5607, 1.9246, 2.3205, 2.7526, 3.1764, 3.5669, 3.8938, 4.1319, 4.2907, 4.3989, 4.5020, 4.5772, 4.6375, 4.7017, 4.7354, 4.7297, 4.7072, 4.6585, 4.6272, 4.6272, 4.6427, 4.6532, 4.6638, 4.6798, 4.6907, 4.7184, 4.7411, 4.7241, 4.6907, 4.6221, 4.5436, 4.4621, 4.3827, 4.3055, 4.2200, 4.1102, 3.9496, 3.6997, 3.3463, 2.9629, 2.5823, 2.2434, 1.9852, 1.8006, 1.6917, 1.6702, 1.7085, 1.8043, 1.9557, 2.1270, 2.3287, 2.5446, 2.7511, 2.9580, 3.1573, 3.3448, 3.5235, 3.6916, 3.8021, 3.8889, 3.9574, 3.9948, 4.0509, 4.0951, 4.1382, 4.1964, 4.2338, 4.2583, 4.2654, 4.2372, 4.2338, 4.2234, 4.1669, 4.1164, 4.0058, 3.8671, 3.7475, 3.5616, 3.3448, 3.1079, 2.8051, 2.4796, 2.1300, 1.7628, 1.4085, 1.0742, 0.7732, 0.5129, 0.3019, 0.1601, 0.0619, 0.0000) ),
    ( 'ASN', (1.8652, 1.8992, 1.9308, 1.9617, 1.9968, 2.0327, 2.0644, 2.0908, 2.1074, 2.1039, 2.1057, 2.1074, 2.1161, 2.1421, 2.1633, 2.1900, 2.2226, 2.2556, 2.2849, 2.2960, 2.2904, 2.2357, 2.1670, 2.1155, 2.0671, 2.0539, 2.0413, 2.0010, 1.9449, 1.8705, 1.8042, 1.7429, 1.6806, 1.6070, 1.4969, 1.3638, 1.2079, 1.0453, 0.8951, 0.7647, 0.6669, 0.5951, 0.5437, 0.4979, 0.4424, 0.3759, 0.2963, 0.2156, 0.1434, 0.0847, 0.0463, 0.0289, 0.0224, 0.0228, 0.0254, 0.0243, 0.0239, 0.0229, 0.0180, 0.0108, 0.0037, 0.0000, 0.0061, 0.0223, 0.0532, 0.1045, 0.1784, 0.2754, 0.3890, 0.5076, 0.6206, 0.7210, 0.8014, 0.8548, 0.8855, 0.9033, 0.9176, 0.9426, 0.9645, 0.9722, 0.9674, 0.9473, 0.9268, 0.9150, 0.9094, 0.9085, 0.9059, 0.9038, 0.8970, 0.8860, 0.8796, 0.8750, 0.8779, 0.9001, 0.9298, 0.9747, 1.0332, 1.0895, 1.1617, 1.2296, 1.3079, 1.3934, 1.4638, 1.5327, 1.5726, 1.6049, 1.6395, 1.6671, 1.7027, 1.7238, 1.7369, 1.7425, 1.7321, 1.7345, 1.7273, 1.7116, 1.7043, 1.6775, 1.6515, 1.6435, 1.6201, 1.5935, 1.5702, 1.5511, 1.5551, 1.5804, 1.6081, 1.6240, 1.6273, 1.6176, 1.5969, 1.5465, 1.4925, 1.4321, 1.3759, 1.3555, 1.3350, 1.3299, 1.3420, 1.3583, 1.3970, 1.4339, 1.4703, 1.5038) ),
    ( 'PRO', (5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8496, 5.8373, 5.8251, 5.8011, 5.7835, 5.5284, 4.8307, 3.5261, 1.8391, 1.0461, 0.5393, 0.1994, 0.1357, 0.1917, 0.3441, 0.6066, 0.8420, 1.0928, 1.3744, 1.6810, 1.9708, 2.2335, 2.4555, 2.6204, 2.7560, 2.8462, 2.8980, 2.9260, 2.9432, 2.9462, 2.9476, 2.9164, 2.8419, 2.7390, 2.5997, 2.4462, 2.2671, 2.0449, 1.7864, 1.4841, 1.1683, 0.8662, 0.5727, 0.2595, 0.0740, 0.0000, 0.0515, 0.3760, 0.8746, 1.6454, 3.2312, 4.6469, 5.3937, 5.7325, 5.7952, 5.8190, 5.8373, 5.8496, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559, 5.8559) ),
    ( 'GLN', (0.1026, 0.1896, 0.3126, 0.4793, 0.7136, 0.9800, 1.2672, 1.5623, 1.8711, 2.1760, 2.4671, 2.7459, 2.9799, 3.1773, 3.3482, 3.4975, 3.6103, 3.7252, 3.8305, 3.8991, 3.9948, 4.0641, 4.0901, 4.0954, 4.0607, 3.9980, 3.9511, 3.9034, 3.8332, 3.7375, 3.5931, 3.3979, 3.1830, 2.9724, 2.7865, 2.6278, 2.4909, 2.3537, 2.2052, 2.0581, 1.9023, 1.7514, 1.6038, 1.4587, 1.3336, 1.2347, 1.1784, 1.1743, 1.2224, 1.3330, 1.5001, 1.7217, 2.0050, 2.3289, 2.6809, 3.0425, 3.3839, 3.6976, 4.0305, 4.3604, 4.6157, 4.8081, 4.8751, 4.9223, 4.9890, 5.0376, 5.0934, 5.1226, 5.1475, 5.1887, 5.1887, 5.1628, 5.1176, 5.0745, 5.0792, 5.1176, 5.1731, 5.2263, 5.2428, 5.2263, 5.1576, 5.0792, 5.0064, 4.9468, 4.8906, 4.7939, 4.6397, 4.4075, 4.0814, 3.7059, 3.2603, 2.8068, 2.4057, 2.0588, 1.7960, 1.6140, 1.5091, 1.4763, 1.5092, 1.5987, 1.7245, 1.8849, 2.0448, 2.1967, 2.3347, 2.4438, 2.5542, 2.6596, 2.7823, 2.9210, 3.0769, 3.2618, 3.4548, 3.6661, 3.8747, 4.0489, 4.1858, 4.2703, 4.3155, 4.3310, 4.3310, 4.3089, 4.2557, 4.1897, 4.0954, 3.9964, 3.8904, 3.7155, 3.5043, 3.2726, 3.0039, 2.7432, 2.4751, 2.1831, 1.8757, 1.5609, 1.2522, 0.9504, 0.6805, 0.4460, 0.2556, 0.1327, 0.0493, 0.0000) ),
    ( 'ARG', (0.0524, 0.1044, 0.1722, 0.2741, 0.4227, 0.5993, 0.8121, 1.0515, 1.3130, 1.5948, 1.8783, 2.1430, 2.3837, 2.5910, 2.7752, 2.9446, 3.0971, 3.2257, 3.3147, 3.3884, 3.4491, 3.5096, 3.5843, 3.6483, 3.7034, 3.7406, 3.7432, 3.7252, 3.6942, 3.6281, 3.5509, 3.4459, 3.3306, 3.2397, 3.1425, 3.0614, 2.9501, 2.7921, 2.6277, 2.4458, 2.2695, 2.1114, 1.9650, 1.8508, 1.7860, 1.7690, 1.8105, 1.9087, 2.0436, 2.2293, 2.4518, 2.6998, 3.0011, 3.3158, 3.6562, 4.0039, 4.2953, 4.5432, 4.7083, 4.8762, 5.0521, 5.2191, 5.3746, 5.4576, 5.4918, 5.5018, 5.4968, 5.4868, 5.4721, 5.4672, 5.4672, 5.4819, 5.5170, 5.5428, 5.5587, 5.5802, 5.5747, 5.5428, 5.5119, 5.4480, 5.3924, 5.3484, 5.2938, 5.2617, 5.2041, 5.0946, 4.9603, 4.7454, 4.4870, 4.1657, 3.8030, 3.4261, 3.0680, 2.7742, 2.5328, 2.3669, 2.2779, 2.2468, 2.2803, 2.3637, 2.4749, 2.6170, 2.7808, 2.9388, 3.0903, 3.2231, 3.3175, 3.4186, 3.5325, 3.6257, 3.7100, 3.7660, 3.7839, 3.8085, 3.8481, 3.8715, 3.9147, 3.9481, 3.9664, 3.9599, 3.9209, 3.8695, 3.7794, 3.7034, 3.6037, 3.4953, 3.3926, 3.2651, 3.1149, 2.9307, 2.7090, 2.4703, 2.2238, 1.9657, 1.6910, 1.4069, 1.1154, 0.8376, 0.5937, 0.3852, 0.2234, 0.1175, 0.0449, 0.0000) ),
    ( 'TRP', (3.8401, 3.8380, 3.8249, 3.8310, 3.8424, 3.8598, 3.8895, 3.8955, 3.8716, 3.8197, 3.7544, 3.6833, 3.6078, 3.5125, 3.4030, 3.2779, 3.1265, 2.9612, 2.7419, 2.4942, 2.2422, 1.9984, 1.7811, 1.6057, 1.4587, 1.3176, 1.1919, 1.0503, 0.9200, 0.7943, 0.6764, 0.5821, 0.4994, 0.4580, 0.4478, 0.4781, 0.5487, 0.6466, 0.7873, 0.9556, 1.1564, 1.3896, 1.6294, 1.8698, 2.0808, 2.2457, 2.3629, 2.4518, 2.5431, 2.6112, 2.6541, 2.6611, 2.6281, 2.5750, 2.4957, 2.3827, 2.2539, 2.1190, 2.0030, 1.8940, 1.8033, 1.7215, 1.6548, 1.6212, 1.5984, 1.6008, 1.6064, 1.5942, 1.5990, 1.6045, 1.6175, 1.6606, 1.6829, 1.7078, 1.7346, 1.7473, 1.7907, 1.8176, 1.8531, 1.9056, 1.9341, 1.9894, 2.0326, 2.0516, 2.0778, 2.0838, 2.0818, 2.0818, 2.0680, 2.0555, 2.0345, 2.0075, 1.9636, 1.8981, 1.8290, 1.7319, 1.6357, 1.5137, 1.3715, 1.2124, 1.0143, 0.8153, 0.6096, 0.4202, 0.2587, 0.1293, 0.0436, 0.0000, 0.0102, 0.0598, 0.1369, 0.2267, 0.3045, 0.3770, 0.4483, 0.5351, 0.6447, 0.7773, 0.9444, 1.1424, 1.3730, 1.6312, 1.9047, 2.1738, 2.4433, 2.7025, 2.9261, 3.1071, 3.2300, 3.3455, 3.4255, 3.5043, 3.5591, 3.5943, 3.6449, 3.6688, 3.7132, 3.7336, 3.7544, 3.8030, 3.8253, 3.8493, 3.8594, 3.8401) ),
    ( 'TYR', (4.3327, 4.3412, 4.3926, 4.4065, 4.4220, 4.4571, 4.4315, 4.4378, 4.4571, 4.4378, 4.4474, 4.4571, 4.4604, 4.5036, 4.5488, 4.5524, 4.5139, 4.3972, 4.2470, 4.0936, 3.9392, 3.8346, 3.7367, 3.6404, 3.5486, 3.3994, 3.2109, 3.0048, 2.7927, 2.5843, 2.3895, 2.1929, 1.9593, 1.6597, 1.3305, 1.0262, 0.7763, 0.6176, 0.5470, 0.5504, 0.6193, 0.7233, 0.8422, 0.9614, 1.0757, 1.1943, 1.3201, 1.4583, 1.6100, 1.7685, 1.9288, 2.0876, 2.2307, 2.3465, 2.4435, 2.5271, 2.6019, 2.6708, 2.7175, 2.7418, 2.7559, 2.7781, 2.8069, 2.8574, 2.9211, 2.9845, 3.0684, 3.1449, 3.2183, 3.2956, 3.3617, 3.4245, 3.4755, 3.5001, 3.4976, 3.4816, 3.4538, 3.4467, 3.4514, 3.4467, 3.4314, 3.3804, 3.3182, 3.2636, 3.2137, 3.1608, 3.0627, 2.9253, 2.7541, 2.5466, 2.3477, 2.1581, 1.9644, 1.7861, 1.5984, 1.3948, 1.1943, 0.9838, 0.7778, 0.5883, 0.4106, 0.2612, 0.1407, 0.0479, 0.0000, 0.0049, 0.0692, 0.2038, 0.3986, 0.6207, 0.8515, 1.0574, 1.2281, 1.3951, 1.5501, 1.7090, 1.8604, 2.0122, 2.1607, 2.2910, 2.4306, 2.5572, 2.6928, 2.8271, 2.9281, 3.0155, 3.0814, 3.1484, 3.2354, 3.3007, 3.3639, 3.3960, 3.4153, 3.4467, 3.4646, 3.5063, 3.5526, 3.6149, 3.6904, 3.7656, 3.8294, 3.8708, 3.9111, 3.9242) ),
  ]

def make_restraints(atmsel, restraints, num_selected):
    from modeller import forms, physical, features
    for (res, values) in chi2:
        arr = True
        for a in atmsel.find_chi2_dihedrals(res, num_selected):
            r = forms.spline(physical.chi2_dihedral,
                             features.dihedral(*a), open=False, low=-3.11978,
                             high=3.16341, delta=0.04363, lowderiv=0,
                             highderiv=0, values=values, use_array=arr)
            arr = restraints.add(r)