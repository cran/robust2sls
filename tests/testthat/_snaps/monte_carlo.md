# generate_param() works correctly

    $structural
    $structural$mean
         [,1]
    [1,] 0.00
    [2,] 0.00
    [3,] 0.00
    [4,] 0.74
    [5,] 0.13
    [6,] 0.66
    [7,] 0.71
    [8,] 0.46
    
    $structural$cov
           [,1]   [,2]   [,3]    [,4]    [,5]    [,6]    [,7]    [,8]
    [1,] 1.0000 0.2636 0.1044  0.0000  0.0000  0.0000  0.0000  0.0000
    [2,] 0.2636 0.8765 0.4356  0.0000  0.0000  0.0000  0.0000  0.0000
    [3,] 0.1044 0.4356 1.1933  0.0000  0.0000  0.0000  0.0000  0.0000
    [4,] 0.0000 0.0000 0.0000  2.0718 -1.3645 -0.6586 -0.0775  1.4605
    [5,] 0.0000 0.0000 0.0000 -1.3645  3.0173  0.3181 -0.2227  0.4211
    [6,] 0.0000 0.0000 0.0000 -0.6586  0.3181  1.1973 -1.1569 -0.2533
    [7,] 0.0000 0.0000 0.0000 -0.0775 -0.2227 -1.1569  2.1569 -0.3381
    [8,] 0.0000 0.0000 0.0000  1.4605  0.4211 -0.2533 -0.3381  2.1077
    
    
    $params
    $params$beta
         [,1]
    [1,] 0.01
    [2,] 0.21
    [3,] 0.91
    [4,] 0.61
    [5,] 0.38
    
    $params$sigma
    [1] 1
    
    $params$Pi
         [,1] [,2] [,3]  [,4]  [,5]
    [1,]    1    0    0 -0.13 -0.93
    [2,]    0    1    0  0.95 -0.14
    [3,]    0    0    1  0.92  0.78
    [4,]    0    0    0  0.14  0.19
    [5,]    0    0    0  0.19  1.00
    [6,]    0    0    0 -0.31 -0.20
    
    $params$Omega2
         [,1]
    [1,] 0.28
    [2,] 0.04
    
    $params$Omega
         [,1]
    [1,] 0.00
    [2,] 0.00
    [3,] 0.00
    [4,] 0.28
    [5,] 0.04
    
    $params$Sigma2_half
         [,1] [,2]
    [1,] 0.91 0.22
    [2,] 0.22 1.07
    
    $params$Sigma_half
         [,1] [,2] [,3] [,4] [,5]
    [1,]    0    0    0 0.00 0.00
    [2,]    0    0    0 0.00 0.00
    [3,]    0    0    0 0.00 0.00
    [4,]    0    0    0 0.91 0.22
    [5,]    0    0    0 0.22 1.07
    
    $params$mean_z
         [,1]
    [1,] 0.74
    [2,] 0.13
    [3,] 0.66
    [4,] 0.71
    [5,] 0.46
    
    $params$cov_z
            [,1]    [,2]    [,3]    [,4]    [,5]
    [1,]  2.0718 -1.3645 -0.6586 -0.0775  1.4605
    [2,] -1.3645  3.0173  0.3181 -0.2227  0.4211
    [3,] -0.6586  0.3181  1.1973 -1.1569 -0.2533
    [4,] -0.0775 -0.2227 -1.1569  2.1569 -0.3381
    [5,]  1.4605  0.4211 -0.2533 -0.3381  2.1077
    
    $params$Ezz
            [,1]    [,2]    [,3]    [,4]    [,5]
    [1,]  2.6194 -1.2683 -0.1702  0.4479  1.8009
    [2,] -1.2683  3.0342  0.4039 -0.1304  0.4809
    [3,] -0.1702  0.4039  1.6329 -0.6883  0.0503
    [4,]  0.4479 -0.1304 -0.6883  2.6610 -0.0115
    [5,]  1.8009  0.4809  0.0503 -0.0115  2.3193
    
    $params$Mzz
         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
    [1,] 1.00  0.7400  0.1300  0.6600  0.7100  0.4600
    [2,] 0.74  2.6194 -1.2683 -0.1702  0.4479  1.8009
    [3,] 0.13 -1.2683  3.0342  0.4039 -0.1304  0.4809
    [4,] 0.66 -0.1702  0.4039  1.6329 -0.6883  0.0503
    [5,] 0.71  0.4479 -0.1304 -0.6883  2.6610 -0.0115
    [6,] 0.46  1.8009  0.4809  0.0503 -0.0115  2.3193
    
    $params$Mxx_tilde_inv
              [,1]       [,2]       [,3]       [,4]       [,5]
    [1,]  27.77782   56.56671   58.01835  -94.82079  13.111081
    [2,]  56.56671  125.03392  126.69627 -205.95344  28.997432
    [3,]  58.01835  126.69627  129.60090 -209.90574  28.929967
    [4,] -94.82079 -205.95344 -209.90574  341.61103 -47.562665
    [5,]  13.11108   28.99743   28.92997  -47.56266   7.201359
    
    
    $setting
    $setting$call
    p <- generate_param(3, 2, 3)
    
    $setting$intercept
    [1] TRUE
    
    $setting$formula
    y ~ x1 + x2 + x3 + x4 + x5 | x1 + x2 + x3 + z4 + z5 + z6
    NULL
    
    $setting$dx1
    [1] 3
    
    $setting$dx2
    [1] 2
    
    $setting$dz2
    [1] 3
    
    
    $names
    $names$x1
    [1] "x1" "x2" "x3"
    
    $names$x2
    [1] "x4" "x5"
    
    $names$x
    [1] "x1" "x2" "x3" "x4" "x5"
    
    $names$z2
    [1] "z4" "z5" "z6"
    
    $names$z
    [1] "x1" "x2" "x3" "z4" "z5" "z6"
    
    $names$r
    [1] "r1" "r2" "r3" "r4" "r5"
    
    $names$u
    [1] "u"
    
    

# generate_data() works correctly

    $data
                y x1          x2          x3           x4          x5           u
    1  -4.3814575  1  2.79200008 -1.68765944 -0.526468959 -3.11393385 -1.93756650
    2   2.0621646  1  0.37752824  1.59231751  2.660881783  0.75772310 -1.38719797
    3   1.3716163  1  1.74991989  0.84521409  0.625886342 -0.61523607  0.07698732
    4  -0.7751130  1  1.33428895  0.01704238  0.388005503  0.68004249 -1.57592170
    5   0.6584386  1  1.51733747 -0.63980232  1.588996280 -1.25433177  0.41937613
    6   2.5424517  1 -0.07879376  0.29486621  2.186068450  3.14610043 -0.24834974
    7  -0.7030792  1  2.32970741 -2.22669475  0.969197495  0.54087056  0.02723322
    8   2.2759349  1  0.17230683  0.46517553  2.329883347  4.16097962 -1.19596036
    9  -3.6955411  1  2.35786773 -3.40430129 -0.126033819  0.45684250 -1.19949862
    10  1.6110858  1  0.87401401  1.03115136  0.387177886  0.24160512  0.15120670
    11 -4.1004764  1  2.51528070 -2.37056372 -0.878377240 -4.74074826 -0.14417787
    12 -4.8279712  1  3.65856766 -3.20495274 -0.025675334 -3.93083337 -1.18038480
    13  5.8137778  1 -0.72301279  3.02658744  2.588952233  2.89125369  0.52347868
    14 -2.7031464  1  0.07389121 -0.68982725 -1.234433211 -3.45539008 -0.03486831
    15  0.6086462  1 -0.08082488 -1.03650453 -0.137900929 -0.08988544  1.67711461
    16  0.2097123  1  1.74971900 -0.02333581  0.262346571  0.83421067 -0.62352457
    17  2.1412310  1 -0.05239055  0.53650174  1.703611868  1.98115764 -0.13802667
    18  7.6953314  1 -2.12061007  5.02837056  2.339131821  2.49184483  1.18107082
    19  1.5805289  1 -2.96220184  2.11963789 -0.210860088  2.25655015 -0.46514364
    20 -2.8657518  1  2.75691319 -2.11453980 -0.520453596 -3.50128664  0.11749329
    21  4.3249571  1  0.76022625  2.27240844  2.613172427  1.60511439 -0.11656074
    22  2.2777726  1 -2.07645729  2.11468933  0.081329510  3.48458290 -0.59429112
    23  2.7199610  1  1.22051545  1.44505796  3.098353085 -0.48499568 -0.56704702
    24 -0.2870096  1  2.50944882 -1.34225958  0.002810355 -0.62041272  0.63150489
    25 -2.1958206  1  3.39169423 -3.19557432  0.472866779 -2.53278041  0.66390409
    26  4.0522315  1  0.44499334  1.54649418  2.098584743  0.46630401  1.08414095
    27  1.6231637  1  0.77415842  0.83928651  2.347928394 -0.61404058 -0.51206114
    28 -0.5813740  1 -1.49711893  2.49286068 -1.235086095 -0.74401636 -1.50935350
    29  1.2594808  1 -0.05028370 -1.76036404  0.323469067  3.43070550  1.36098744
    30  2.5025444  1  0.09877949  1.91096065  1.575087691  1.91713665 -0.95648895
       x1          x2          x3          z4          z5          z6 r1 r2 r3
    1   1  2.79200008 -1.68765944  0.39685269  0.32924413  1.69526291  0  0  0
    2   1  0.37752824  1.59231751  0.83062951 -0.11838824  0.91335925  0  0  0
    3   1  1.74991989  0.84521409  0.54208177  0.51815124  2.29374880  0  0  0
    4   1  1.33428895  0.01704238 -0.84570974  3.00329350  0.96146249  0  0  0
    5   1  1.51733747 -0.63980232  1.03903520 -0.49421663  0.86469942  0  0  0
    6   1 -0.07879376  0.29486621 -1.44504138  3.28686838 -1.20990445  0  0  0
    7   1  2.32970741 -2.22669475 -0.42949125  2.44781374  0.83706218  0  0  0
    8   1  0.17230683  0.46517553 -0.14777493  2.73077381  0.09000197  0  0  0
    9   1  2.35786773 -3.40430129 -2.39297556  4.31818488 -0.76647916  0  0  0
    10  1  0.87401401  1.03115136 -0.50551960  0.93517383  0.49171853  0  0  0
    11  1  2.51528070 -2.37056372  0.86792371 -0.13224850  1.02534043  0  0  0
    12  1  3.65856766 -3.20495274 -0.83800590  1.51495634  1.37241757  0  0  0
    13  1 -0.72301279  3.02658744  1.35626776  0.29962133  0.84282138  0  0  0
    14  1  0.07389121 -0.68982725  0.73354851 -0.02235288 -1.36516439  0  0  0
    15  1 -0.08082488 -1.03650453  0.24519669  1.15463472 -1.71904774  0  0  0
    16  1  1.74971900 -0.02333581  0.47966090  1.38767954  1.91550764  0  0  0
    17  1 -0.05239055  0.53650174  0.11942244  2.31985091 -0.16948468  0  0  0
    18  1 -2.12061007  5.02837056  2.57869064 -0.79026783  0.82499578  0  0  0
    19  1 -2.96220184  2.11963789  3.32906525  0.05796597 -1.68915909  0  0  0
    20  1  2.75691319 -2.11453980  1.02517328 -0.64573396  1.52155894  0  0  0
    21  1  0.76022625  2.27240844 -1.04181951  1.28147103  1.01325893  0  0  0
    22  1 -2.07645729  2.11468933  1.32795224  1.32229598 -1.37708335  0  0  0
    23  1  1.22051545  1.44505796  1.75926693 -0.42425775  2.48384335  0  0  0
    24  1  2.50944882 -1.34225958  1.17301350  0.91343288  2.36437712  0  0  0
    25  1  3.39169423 -3.19557432 -0.06164772 -0.46907700  0.83861174  0  0  0
    26  1  0.44499334  1.54649418  0.05361181  0.61376268  0.58325237  0  0  0
    27  1  0.77415842  0.83928651  0.90379787 -0.45357691  0.83319461  0  0  0
    28  1 -1.49711893  2.49286068  1.37856659 -0.37341714 -0.88630618  0  0  0
    29  1 -0.05028370 -1.76036404 -1.65057767  4.50305169 -2.38719734  0  0  0
    30  1  0.09877949  1.91096065  0.43057923  1.04630635  0.77509271  0  0  0
                r4          r5
    1  -1.08880661 -0.54227304
    2   1.15664485  0.64180988
    3  -1.14741257 -0.26191129
    4  -0.91942101 -0.86676626
    5   1.08263687  0.85688111
    6   1.32237602  0.58078225
    7   0.78906806  1.33505353
    8   1.39797245  2.06756948
    9   0.17290193  0.35548363
    10 -1.21627222 -0.25111238
    11 -0.73550186 -1.43715819
    12 -0.16782952 -1.07002246
    13  0.63582302  0.97054576
    14 -1.06163950 -2.36703427
    15  0.23585395  0.09224098
    16 -0.98542180  0.93165967
    17  0.87977023  0.10891350
    18 -0.09750664 -0.33185385
    19 -0.21755719  0.09020416
    20 -0.61329604  0.21928508
    21  0.14682723 -0.01180605
    22 -0.62559564  0.62439763
    23  1.34371312  0.07549698
    24 -0.62110444  1.04444256
    25  0.67841059  2.01311716
    26  0.43975409 -0.25496084
    27  1.18827255  0.21819236
    28 -2.37305992 -2.03381609
    29  0.75624338  1.05986833
    30 -0.16563718  0.39731861
    
    $beta
         [,1]
    [1,] 0.01
    [2,] 0.21
    [3,] 0.91
    [4,] 0.61
    [5,] 0.38
    
    $Pi
         [,1] [,2] [,3]  [,4]  [,5]
    [1,]    1    0    0 -0.13 -0.93
    [2,]    0    1    0  0.95 -0.14
    [3,]    0    0    1  0.92  0.78
    [4,]    0    0    0  0.14  0.19
    [5,]    0    0    0  0.19  1.00
    [6,]    0    0    0 -0.31 -0.20
    

# mc_grid() works correctly

         M    n iterations sign_level initial_est split mean_gauge        avar
    1  100  100          0       0.01 robustified   0.5    0.00820 0.007125496
    2  100 1000          0       0.01 robustified   0.5    0.00969 0.007125496
    3  100  100          0       0.05 robustified   0.5    0.04870 0.021256489
    4  100 1000          0       0.05 robustified   0.5    0.05054 0.021256489
    5  100  100          0       0.01   saturated   0.3    0.03350 0.009239404
    6  100 1000          0       0.01   saturated   0.3    0.01182 0.009239404
    7  100  100          0       0.05   saturated   0.3    0.10020 0.041251545
    8  100 1000          0       0.05   saturated   0.3    0.05525 0.041251545
    9  100  100          0       0.01   saturated   0.4    0.02570 0.007587914
    10 100 1000          0       0.01   saturated   0.4    0.01107 0.007587914
    11 100  100          0       0.05   saturated   0.4    0.08660 0.025630407
    12 100 1000          0       0.05   saturated   0.4    0.05379 0.025630407
    13 100  100          0       0.01   saturated   0.5    0.02290 0.007125496
    14 100 1000          0       0.01   saturated   0.5    0.01070 0.007125496
    15 100  100          0       0.05   saturated   0.5    0.08060 0.021256489
    16 100 1000          0       0.05   saturated   0.5    0.05283 0.021256489
       mean_avar_est    var_gauge var_ratio var_ratio2 mean_prop_t mean_prop_p
    1    0.007125496 5.531313e-05 0.7762706  0.7762706   0.5923282   0.6398972
    2    0.007125496 6.438283e-06 0.9035557  0.9035557   0.7679736   0.5073626
    3    0.021256489 2.195051e-04 1.0326496  1.0326496   0.8024909   0.4971735
    4    0.021256489 2.057414e-05 0.9678993  0.9678993   0.7678166   0.5184822
    5    0.009239404 5.219697e-04 5.6493872  5.6493872   2.5696559   0.2733855
    6    0.009239404 1.394707e-05 1.5095206  1.5095206   1.0856553   0.4045535
    7    0.041251545 1.642384e-03 3.9813875  3.9813875   2.5011723   0.2053616
    8    0.041251545 6.313889e-05 1.5305824  1.5305824   1.1350302   0.4161983
    9    0.007587914 2.873838e-04 3.7873894  3.7873894   1.9860259   0.2851764
    10   0.007587914 1.065162e-05 1.4037609  1.4037609   0.9475002   0.4615559
    11   0.025630407 1.044889e-03 4.0767549  4.0767549   2.4235611   0.2193477
    12   0.025630407 3.687465e-05 1.4387070  1.4387070   1.1436702   0.3725455
    13   0.007125496 3.218081e-04 4.5162900  4.5162900   1.7651382   0.3742689
    14   0.007125496 9.707071e-06 1.3623010  1.3623010   0.9140759   0.4653316
    15   0.021256489 9.349899e-04 4.3986093  4.3986093   2.3457425   0.2006471
    16   0.021256489 3.743545e-05 1.7611307  1.7611307   1.0996695   0.4196612
       prop_size_001 prop_size_005 prop_size_010 mean_count_p count_size_001
    1           0.00          0.04          0.04    0.9043513           0.00
    2           0.01          0.02          0.05    0.6140676           0.00
    3           0.02          0.06          0.06    0.7427102           0.00
    4           0.02          0.04          0.07    0.6713238           0.00
    5           0.45          0.57          0.57    0.3181477           0.28
    6           0.06          0.18          0.26    0.4487820           0.06
    7           0.38          0.57          0.57    0.2380981           0.30
    8           0.07          0.21          0.27    0.4571812           0.05
    9           0.23          0.45          0.45    0.3637978           0.13
    10          0.03          0.09          0.14    0.5483525           0.03
    11          0.37          0.47          0.62    0.3009331           0.16
    12          0.07          0.16          0.18    0.4994845           0.00
    13          0.21          0.33          0.33    0.4681521           0.13
    14          0.05          0.08          0.11    0.5752778           0.02
    15          0.40          0.54          0.54    0.3309187           0.12
    16          0.10          0.15          0.21    0.5626641           0.01
       count_size_005 count_size_010
    1            0.00           0.04
    2            0.01           0.02
    3            0.01           0.01
    4            0.01           0.02
    5            0.45           0.57
    6            0.12           0.18
    7            0.46           0.46
    8            0.14           0.24
    9            0.23           0.45
    10           0.07           0.09
    11           0.37           0.38
    12           0.06           0.12
    13           0.21           0.33
    14           0.05           0.08
    15           0.30           0.30
    16           0.01           0.10
