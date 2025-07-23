# Drip rate estimation using trace element concentrations in speleothem
# The script processes data from an Excel file to estimate drip rates 
# using trace element concentrations found in speleothems. 
#

# Import Libraries
import os, sys, io, openpyxl, copy, pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import matplotlib.ticker as ticker
from PIL import Image

from drip_rate_util import *

from params import OUTLIER_WINDOW_SIZE as WS
from params import OUTLIER_ZERO_TOL as ZERO_TOL
from params import OUTLIER_PDF_TOL as PDF_TOL


# Initialization of settings
bErrorCheck = False
bOutlierProcessing = True

# Load Excel sheets
xlName = 'Drip_rate.xlsx'

wb = openpyxl.load_workbook(xlName)

nameSheets = ['1.Depth_Age', '2.Trace_Elems','3.Isotopes','0.Settings']
wsheets = []
bError = False
for name in nameSheets:
    try:
        wsheets.append(wb[name])
    except:
        print("Error: {} sheet is missing".format(name))
if bErrorCheck and bError:
    exit()

# Extract settings from the Excel sheet
StationName = wsheets[3].cell(2,2).value
StationNameShort = wsheets[3].cell(3,2).value

# Starting row in the Excel sheet
baseRow = 21
# Setting for proxy record calculation
# Proxy record calculation takes a long time (up to 1 hour)
# print(wsheets[3].cell(baseRow,2).value)
bProxyRecord = True if wsheets[3].cell(baseRow,2).value in ['Yes','True',1] else False
# print(bProxyRecord)
# bCalcIsotope = True if wsheets[3].cell(baseRow+1,2).value in ['Yes','True',1] else False
# bGenPlot = True if wsheets[3].cell(baseRow+1,2).value in ['Yes','True',1] else False

# Extract Plot configuration
plotConfs =[]
baseRow = 25
plotConf = {}
val = wsheets[3].cell(baseRow,3).value
plotConf['Width'] = val/25.4 if val is not None and type(val) in [float, int] else 14.0
val = wsheets[3].cell(baseRow+1,3).value
plotConf['Height'] = val/25.4 if val is not None and type(val) in [float, int] else 5.0
val = wsheets[3].cell(baseRow+2,3).value
plotConf['Xlabel'] = val if val and type(val) in [str] else 'X-label'
val = wsheets[3].cell(baseRow+3,3).value
plotConf['Xmin'] = float(val) if val is not None and type(val) in [float, int] else -9500.0
val = wsheets[3].cell(baseRow+4,3).value
plotConf['Xmax'] = float(val) if val is not None and type(val) in [float, int] else 100.0
val = wsheets[3].cell(baseRow+5,3).value
plotConf['XtickMajor'] = float(val) if val is not None and type(val) in [float, int] else 200.0
val = wsheets[3].cell(baseRow+6,3).value
plotConf['XtickMinor'] = float(val) if val and type(val) in [float, int] else 100
val = wsheets[3].cell(baseRow+7,3).value
plotConf['Ylabel'] = val if val and type(val) in [str] else 'Y-label'
val = wsheets[3].cell(baseRow+8,3).value
plotConf['Ymin'] = float(val) if val is not None and type(val) in [float, int] else 0.0
val = wsheets[3].cell(baseRow+9,3).value
plotConf['Ymax'] = float(val) if val is not None and type(val) in [float, int] else 35.0
val = wsheets[3].cell(baseRow+10,3).value
plotConf['YtickMajor'] = float(val) if val is not None and type(val) in [float, int] else 1.0
val = wsheets[3].cell(baseRow+11,3).value
plotConf['YtickMinor'] = float(val) if val and type(val) in [float, int] else 100
val = wsheets[3].cell(baseRow+12,3).value
plotConf['Plot1Clr'] = val if val is not None and type(val) in [str] else 'black'
val = wsheets[3].cell(baseRow+13,3).value
plotConf['Plot1Type'] = val if val and type(val) in [str] else '-'
val = wsheets[3].cell(baseRow+14,3).value
plotConf['Plot1Label'] = val if val and type(val) in [str] else ''
val = wsheets[3].cell(baseRow+15,3).value
plotConf['Plot2Clr'] = val if val is not None and type(val) in [str] else 'black'
val = wsheets[3].cell(baseRow+16,3).value
plotConf['Plot2Type'] = val if val and type(val) in [str] else '-'
val = wsheets[3].cell(baseRow+17,3).value
plotConf['Plot2Label'] = val if val and type(val) in [str] else ''
val = wsheets[3].cell(baseRow+18,3).value
plotConf['Plot3Clr'] = val if val is not None and type(val) in [str] else 'black'
val = wsheets[3].cell(baseRow+19,3).value
plotConf['Plot3Type'] = val if val and type(val) in [str] else '-'
val = wsheets[3].cell(baseRow+20,3).value
plotConf['Plot3Label'] = val if val and type(val) in [str] else ''
val = wsheets[3].cell(baseRow+21,3).value
plotConf['CMapLabel'] = val if val and type(val) in [str] else ''
plotConfs.append(plotConf)

baseRow = 50
plotConf = {}
val = wsheets[3].cell(baseRow,3).value
plotConf['Width'] = val/25.4 if val is not None and type(val) in [float, int] else 14.0
val = wsheets[3].cell(baseRow+1,3).value
plotConf['Height'] = val/25.4 if val is not None and type(val) in [float, int] else 5.0
val = wsheets[3].cell(baseRow+2,3).value
plotConf['Xlabel'] = val if val and type(val) in [str] else 'X-label'
val = wsheets[3].cell(baseRow+3,3).value
plotConf['Xmin'] = float(val) if val is not None and type(val) in [float, int] else -9500.0
val = wsheets[3].cell(baseRow+4,3).value
plotConf['Xmax'] = float(val) if val is not None and type(val) in [float, int] else 100.0
val = wsheets[3].cell(baseRow+5,3).value
plotConf['XtickMajor'] = float(val) if val is not None and type(val) in [float, int] else 200.0
val = wsheets[3].cell(baseRow+6,3).value
plotConf['XtickMinor'] = float(val) if val and type(val) in [float, int] else 100
val = wsheets[3].cell(baseRow+7,3).value
plotConf['Y1label'] = val if val and type(val) in [str] else 'Y1-label'
val = wsheets[3].cell(baseRow+8,3).value
plotConf['Y1min'] = float(val) if val is not None and type(val) in [float, int] else 0.0
val = wsheets[3].cell(baseRow+9,3).value
plotConf['Y1max'] = float(val) if val is not None and type(val) in [float, int] else 35.0
val = wsheets[3].cell(baseRow+10,3).value
plotConf['Y1tickMajor'] = float(val) if val is not None and type(val) in [float, int] else 1.0
val = wsheets[3].cell(baseRow+11,3).value
plotConf['Y1tickMinor'] = float(val) if val and type(val) in [float, int] else 100
val = wsheets[3].cell(baseRow+12,3).value
plotConf['Plot1Clr'] = val if val is not None and type(val) in [str] else 'black'
val = wsheets[3].cell(baseRow+13,3).value
plotConf['Plot1Type'] = val if val and type(val) in [str] else '-'
val = wsheets[3].cell(baseRow+14,3).value
plotConf['Plot1Label'] = val if val and type(val) in [str] else ''
val = wsheets[3].cell(baseRow+15,3).value
plotConf['Y2label'] = val if val and type(val) in [str] else 'Y2-label'
val = wsheets[3].cell(baseRow+16,3).value
plotConf['Y2min'] = float(val) if val is not None and type(val) in [float, int] else -20.0
val = wsheets[3].cell(baseRow+17,3).value
plotConf['Y2max'] = float(val) if val is not None and type(val) in [float, int] else 20.0
val = wsheets[3].cell(baseRow+18,3).value
plotConf['Y2tickMajor'] = float(val) if val is not None and type(val) in [float, int] else 1.0
val = wsheets[3].cell(baseRow+19,3).value
plotConf['Y2tickMinor'] = float(val) if val and type(val) in [float, int] else 100
val = wsheets[3].cell(baseRow+20,3).value
plotConf['Plot2Clr'] = val if val is not None and type(val) in [str] else 'black'
val = wsheets[3].cell(baseRow+21,3).value
plotConf['Plot2Type'] = val if val and type(val) in [str] else '-'
val = wsheets[3].cell(baseRow+22,3).value
plotConf['Plot2Label'] = val if val and type(val) in [str] else ''
val = wsheets[3].cell(baseRow+23,3).value
plotConf['Plot2Isotope'] = int(val) if val and type(val) in [int] else 1
plotConfs.append(plotConf)

nameSheets = ['4.Output','5.OutDripRate','6.OutIsotope', '7.PlotDripRate','8.PlotDripIsotope']
outWSs = []

# wb.save(xlName)


# If the proxy record calculation is required,
if bProxyRecord:
    # Check data column settings
    for name in nameSheets:
        if name in wb.sheetnames:
            del(wb[name])
        outWSs.append(wb.create_sheet(name))
    print('Start a long proxy record processing. This will take more than 40 minutes')
    # , ['Depth Output', 'Age Output', 'Error Output']
    # Check Depth Age
    paramNames = [[[['Depth Input', 'Age Input', 'Error Input']]],
                [[['Depth1 Input', 'Trace Element1 Input'], ['Depth1 Output', 'Trace Element1 Output']],
                [['Depth2 Input', 'Trace Element2 Input'], ['Depth2 Output', 'Trace Element2 Output']]],
                [[['Depth(d18O) Input','d18O Input'], ['Depth(d18O) Output', 'd18O Output']], 
                [['Depth(d13C) Input', 'd13C Input'], ['Depth(d13C) Output', 'd13C Output']]]
                ]
    cellNames = [[[['B5', 'B6', 'B7']]],
                [[['B5', 'B6'], ['C5', 'C6']], [['B7', 'B8'], ['C7', 'C8']]],
                [[['B5', 'B6'], ['C5', 'C6']], [['B7', 'B8'], ['C7', 'C8']]]]

    WindowSizeAdd = 'B9'

    data_lst = []
    iSheet = 0
    bError = False
    data_lst.append([])
    for idx1, list1 in enumerate(paramNames[iSheet]):
        data_lst[iSheet].append([])
        # break
        for idx2, row2 in enumerate(list1):
            data_lst[iSheet][idx1].append([])
            # break
            for idx3, parName in enumerate(row2):
                # data_lst[iSheet][idx1][idx2].append([])
                # break
                cellName = cellNames[iSheet][idx1][idx2][idx3]
                row1, col1 = openpyxl.utils.cell.coordinate_to_tuple(cellName)
                it1 = wsheets[iSheet].cell(row1,col1).value
                try:
                    row2, col2 = openpyxl.utils.cell.coordinate_to_tuple(it1)
                    data_lst[iSheet][idx1][idx2].append([row2,col2])
                except:
                    bError = True
                    print('Error: {0} ({1}) = "{2}", but "{2}" is not a correct cell reference.'.format(
                        parName, cellName, it1))
    if bErrorCheck and bError:
        exit()



    # Check Trace element settings
    iSheet = 1
    bError = False
    data_lst.append([])
    for idx1, list1 in enumerate(paramNames[iSheet]):
        data_lst[iSheet].append([])
        # break
        for idx2, row2 in enumerate(list1):
            data_lst[iSheet][idx1].append([])
            # break
            for idx3, parName in enumerate(row2):
                # data_lst[iSheet][idx1][idx2].append([])
                # break
                cellName = cellNames[iSheet][idx1][idx2][idx3]
                row1, col1 = openpyxl.utils.cell.coordinate_to_tuple(cellName)
                it1 = wsheets[iSheet].cell(row1,col1).value
                try:
                    row2, col2 = openpyxl.utils.cell.coordinate_to_tuple(it1)
                    data_lst[iSheet][idx1][idx2].append([row2,col2])
                except:
                    bError = True
                    print('Error: {0} ({1}) = "{2}", but "{2}" is not a correct cell reference.'.format(
                        parName, cellName, it1))
    if bErrorCheck and bError:
        exit()



    # Check Isotpe settings
    iSheet = 2
    bError = False
    data_lst.append([])
    for idx1, list1 in enumerate(paramNames[iSheet]):
        data_lst[iSheet].append([])
        # break
        for idx2, row2 in enumerate(list1):
            data_lst[iSheet][idx1].append([])
            # break
            for idx3, parName in enumerate(row2):
                # data_lst[iSheet][idx1][idx2].append([])
                # break
                cellName = cellNames[iSheet][idx1][idx2][idx3]
                row1, col1 = openpyxl.utils.cell.coordinate_to_tuple(cellName)
                it1 = wsheets[iSheet].cell(row1,col1).value
                try:
                    row2, col2 = openpyxl.utils.cell.coordinate_to_tuple(it1)
                    data_lst[iSheet][idx1][idx2].append([row2,col2])
                except:
                    bError = True
                    print('Error: {0} ({1}) = "{2}", but "{2}" is not a correct cell reference.'.format(
                        parName, cellName, it1))
    if bErrorCheck and bError:
        exit()



    # Convert the data into a list
    # DF is a panda dataframe
    # iSheet = 0
    maxRow = 50001
    dataAll = []
    bError = False
    # iSheet = 0
    for iSheet in range(3):
        dataAll.append([])
        # idx1=0;list1=paramNames[iSheet][idx1]
        for idx1, list1 in enumerate(paramNames[iSheet]):
            dataAll[iSheet].append([])
            # break
            # Do not check output
            # idx2=0;row2 = list1[:-1][idx2]
            # lastIdx = -1 if len(list1)>1 else 1
            for idx2, row2 in enumerate(list1):
                if row2[0][-6:]=='Output': continue
                # sys.stop1
                dataAll[iSheet][idx1].append([])
                # break
                name1 = ['age','trace_elem','isotope'][iSheet]
                dataSet = []
                # idx3=0;parName=row2[idx3]
                for idx3, parName in enumerate(row2):
                    list2 = data_lst[iSheet][idx1][idx2][idx3]
                    header1 = wsheets[iSheet].cell(*list2).value
                    name2 = header1.replace('(',' ').split()[0]
                    # quit2()

                    if type(header1)!= str: 
                        bError = True
                        if not "output" in parName.lower(): 
                            print('Error in "{}" sheet: Header for ({}) = "{}" is not string'.format(nameSheets[iSheet], parName, header1))
                    dataSet.append([list2[0]+1,list2[1]])
                # idx4=0
                tLst = []
                nBlankLine = 0
                for idx4 in range(maxRow):
                    bError1 = False
                    # break
                    tLst1 = []
                    for idx3, parName in enumerate(row2):
                        parName = paramNames[iSheet][idx1][idx2][idx3]
                        # break
                        list2 = [dataSet[idx3][0]+idx4,dataSet[idx3][1]]
                        cell1 = wsheets[iSheet].cell(*list2)
                        it1 = cell1.value
                        cellAddress = cell1.coordinate
                        if it1 is None:
                            tLst1.append(None)
                        elif type(it1)==str:
                            try:
                                it1 = float(it1)
                                nBlankLine = 0
                            except:
                                print('Error in "{}" sheet: Not a value "{}" at {}'.format(nameSheets[iSheet], it1, cellAddress))
                                bError = True
                                bError1 = True
                        tLst1.append(it1)
                    if any(tLst1): 
                        tLst.append(tLst1)
                    else:
                        nBlankLine += 1
                    if nBlankLine>=3: break
                # Make dataAll have five entries
                dataAll[iSheet][idx1][idx2] = [[name1,name2], pd.DataFrame(tLst,columns=row2), None, None, None]

    if bErrorCheck and bError:
        exit()

    # Process outliers and replace them as expected median
    # iSheet = 2
    bError = False
    for iSheet in range(1,3):
        # break
        for idx1, list1 in enumerate(paramNames[iSheet]):
            # idx1=0; list1 = paramNames[iSheet][idx1]
            # Do not check output
            for idx2, cols2 in enumerate(list1):
                if cols2[0][-6:]=='Output': continue
                # idx2=0; cols2 = list1[:-1][idx2]
                # break
                # quit2()
                names, df1, dummy = dataAll[iSheet][idx1][idx2][:3]
                name1, name2 = names

                x = df1[cols2[0]].to_numpy()
                y = df1[cols2[1]].to_numpy()
                y_ = y.copy()
                i = detect_outliers(name1, x, y, winsize=WS,
                                        zero_tol=ZERO_TOL, pdf_tol=PDF_TOL)
                
                # replace the outlier with local median excluding itself
                # _printmsg("replace outliers ...", verb)
                if bOutlierProcessing:
                    kk = int(WS / 2)
                    # j = 0
                    for j in range(len(i)):
                        ii = i[j]
                        med = np.median(y[(ii - kk):(ii + kk)])
                        y_[ii] = med

                # save outlier-removed dataset
                # _printmsg("save output ...", verb)
                X = np.c_[x, y_]
                dataAll[iSheet][idx1][idx2][2]=X
                # Output cell
                # Copy header
                inCells = cellNames[iSheet][idx1][0].copy()
                addressi =[]
                for idx3, row3 in enumerate(inCells):
                    add = openpyxl.utils.cell.coordinate_to_tuple(row3)
                    addressi.append(add)
                    add2 = wsheets[iSheet].cell(*add).value
                    addressi[idx3] = openpyxl.utils.cell.coordinate_to_tuple(add2)

                outCells = cellNames[iSheet][idx1][1].copy()
                addresso = []
                for idx3, row3 in enumerate(outCells):
                    add = openpyxl.utils.cell.coordinate_to_tuple(row3)
                    addresso.append(add)
                    add2 = wsheets[iSheet].cell(*add).value
                    addresso[idx3] = openpyxl.utils.cell.coordinate_to_tuple(add2)
                    wsheets[iSheet].cell(*addresso[idx3]).value = wsheets[iSheet].cell(*addressi[idx3]).value

                # Copy data
                # idx4=0
                for idx4, X1 in enumerate(X):
                    # idx5 = 0
                    for idx5, _ in enumerate(X1):
                        # break
                        wsheets[iSheet].cell(addresso[idx5][0]+1+idx4, addresso[idx5][1]).value = X1[idx5]    

    if bErrorCheck and bError:
        exit()
    # End of 001_replace_outliers.py
    # wb.save(xlName)



    # Construct an age dating model
    dtype_lab = "none"
    # print("Creating BayProX:Data:SampleInfo object ...")
    info = dat.SampleInfo(name="none", datatype=dtype_lab, archive="Speleothem",)


    # load age-depth data
    type1 = dataAll[0][0][0][0][0]
    xname = "Distance (cm)"
    yname = "Age(corr) (yrs BP)"
    ename = "Age error (yrs BP 2-sigma)"
    iSheet = 0
    dating_depth = dataAll[iSheet][0][0][1][paramNames[iSheet][0][0][0]].to_numpy()
    dating_age = dataAll[iSheet][0][0][1][paramNames[iSheet][0][0][1]].to_numpy()
    dating_age_error = dataAll[iSheet][0][0][1][paramNames[iSheet][0][0][2]].to_numpy() / 2.

    # set BayProx:Data:DatingTable
    # print("Creating BayProX:Data:DatingTable object ...")
    # pd.Series standard deviation uses N-1 as a degree of freedom
    DT = dat.DatingTable(depth=dating_depth, age=dating_age, 
        ageerror=pd.Series(dating_age_error),
        datingmethod="U/Th", sampleinfo=info)
    DT.calibration()



    # load proxy-depth data for the first trace element series and isotopic
    # ratio series
    type1 = dataAll[1][0][0][0] # First trace element
    xname = "Depth (from top)"
    x_Ni = dataAll[1][0][0][2][:,0]
    xname = "# Depth_cm"
    x_d18O = dataAll[2][0][0][2][:,0] # First isotopic ratio series

    # combine both depth scales into one unified series
    unified_depth = np.sort(np.r_[x_Ni, x_d18O])
    xmin = max(x_Ni.min(), x_d18O.min())
    xmax = min(x_Ni.max(), x_d18O.max())
    i = (unified_depth >= xmin) * (unified_depth <= xmax)
    unified_depth = unified_depth[i]
    # unified_depth = np.unique(unified_depth)

    # Construct a proxy depth relationship using the unified depth series
    PD1 = dat.ProxyDepth(depth=unified_depth, proxy=np.zeros(len(unified_depth)),
        proxyerror=0., sampleinfo=info)


    # choose the prior remainder variance of the age errors such that the
    # standard deviation and mean of the posteriors are as close to the
    # standard deviation and mean of the age measurements
    # print("optimal prior remainder variance ...")
    # This takes 01:00
    verb =-1
    n = 100
    pr_rem_var = np.linspace(0., 0.20, n)
    pr_rem_residual = np.zeros(n)
    # pb = _progressbar_start(max_value=n, pbar_on=True)
    # i = 0
    for i in range(n):
        pr_rem_residual[i] = preremvar_optimize_residual(pr_rem_var[i],
                                            DT, PD1, calage, max_degree,
                                            method, verb
                                            )
    # outWS.cell(1,1).value='Prior Remainder Variance'
    # outWS.cell(1,2).value='Estimated Residual'
    # for idx1 in range(n):
    #     outWS.cell(2+idx1,1).value=pr_rem_var[idx1]
    #     outWS.cell(2+idx1,2).value=pr_rem_residual[idx1]


    # Set the Prior remainder variance series in Excel sheet
    baseRow = 5
    baseCol = 1
    outWSs[0].cell(baseRow,baseCol).value='Prior Remainder Variance'
    outWSs[0].cell(baseRow,baseCol+1).value='Estimated Residual'
    for idx1 in range(n):
        outWSs[0].cell(baseRow+1+idx1,baseCol).value=pr_rem_var[idx1]
        outWSs[0].cell(baseRow+1+idx1,baseCol+1).value=pr_rem_residual[idx1]


    # Store Unified Depth in Excel sheet
    baseRow = 5
    baseCol = 5
    outWSs[0].cell(baseRow,baseCol).value='Unified depth(cm)'
    for idx1 in range(len(x)):
        outWSs[0].cell(baseRow+1+idx1,baseCol).value=x[idx1]


    # Store Proxy Record
    i_proxy_min = np.argmin(pr_rem_residual)
    # set BayProx:AgeDepth:DWF
    # DWF.output contains all the relevant info
    # print("Estimating DWFs ...")
    # Fix MoTaBaR parameters
    pre_remvar = pr_rem_var[i_proxy_min]

    DT1 = dat.DatingTable(depth=dating_depth, age=dating_age, 
        ageerror=dating_age_error,
        datingmethod="U/Th", sampleinfo=info)
    DT1.calibration()

    pre_remvar *= DT1.ageerror.std()
    # pre_remvar *= DT1.ageerror.std()
    # End of 002_bayorixy.py preremvar


    # Process the data for proxyrecord
    # Multiprocessing is needed
    # set BayProx:Data:SampleInfo
    bError = False
    pLst = []
    iTarget = -1
    # iSheet = 1
    for iSheet in range(1,3):
        # break
        dtype_lab = ['age','concentrations','isotpic ratio'][iSheet]

        # print("Creating BayProX:Data:SampleInfo object ...")
        # idx1 = 0; list1=paramNames[iSheet][idx1]
        for idx1, list1 in enumerate(paramNames[iSheet]):
            # Do not check output
            # idx2=0; cols2 = list1[idx2]
            for idx2, cols2 in enumerate(list1):
                if cols2[0][-6:]=='Output': continue
                iTarget += 1
                # quit3()
                # Start of 002_bayproxy.py get_proxyrecord
                names, df1, np1 = dataAll[iSheet][idx1][idx2][:3]# 
                type1, name2 = names
                if iSheet==1:
                    xname = "Depth (from top)"
                    yname = "%s (ppb)" % name2
                elif iSheet == 2:
                    xname = "Depth (cm)"
                    yname = "%s (per mil)" % name2
                # sys.quit2
                # break

                info = dat.SampleInfo(name=StationName,
                    datatype=dtype_lab, archive="Speleothem",)
                # x = df1[paramNames[iSheet][idx1][idx2][0]].to_numpy()
                # y = df1[paramNames[iSheet][idx1][idx2][1]].to_numpy()
                x = np1[:,0]
                y = np1[:,1]
        
                # remove NaNs (there is only on e Nan for Co and Ni at the 29th entry
                i = ~ np.isnan(y)
                x, y = x[i], y[i]


                # interpolate the proxy measurements on to the unified depth scale grid,
                # but only within the depth range on which the proxy was actually measured
                f_y = interp1d(x, y, kind="linear", bounds_error=False, fill_value=np.nan)
                y_ = f_y(unified_depth)


                # set BayProx:Data:ProxyDepth
                # print("Creating BayProX:Data:ProxyDepth object ...")
                # print("\tfor %s ..." % name)
                PD = dat.ProxyDepth(depth=unified_depth, proxy=y_,
                        proxyerror=0., sampleinfo=info)


                # set BayProx:AgeDepth:DWF
                # DWF.output contains all the relevant info
                # print("Estimating DWFs ...")
                # Fix MoTaBaR parameters
                # This takes a while (more than 5 min?)
                # pre_remvar *= DT.ageerror.std()
                DWF = ad.DWF(DT1)
                # Parallel Start
                tLst1 = [iTarget, PD, calage,max_degree,pre_remvar,method]
                pLst.append(tLst1)
#%%
                # took 9:00
                DWF(PD.depth, calBP_axis=calage,
                    max_degree=max_degree,            # maximum degree of regression curve
                    pre_remvar=pre_remvar,            # influences the variance of results
                    method=method, verbose=0)
                # len(DWF.calbp)
#%%    
                # set BayProx:ProxyRecord:ProxyDistribution
                # print("Estimating proxy distributions ...")
                PDist = prx.ProxyDistributions(DWF)
#%%                
                if iSheet == 1:
                    p25 = np.percentile(PD.proxy, 25)
                    p50 = np.percentile(PD.proxy, 50)
                    p75 = np.percentile(PD.proxy, 75)
                    iqr = p75 - p25
                    limits = [0., p50 + 10. * iqr]
                elif iSheet == 2:
                    sd_mult = 3.
                    limits = PDist.get_limits(DWF, PD, sd_mult)
                # Took 1:00
                PDist.get_pdf(DWF, PD, res=500, limits=limits)
                # Parallel End

                # Obtain PDist, DWF, PD

                cdfmat = get_cdfmat(PDist.pdfmat, PDist.proxyspan)
    
                # set BayProx:ProxyRecord:ProxyEstimate
                # print("Estimating proxy mean and variance ...")
                PEst= prx.ProxyEstimates()
                prx_mean= PEst.get_mean(DWF, PD)
                prx_var= PEst.get_variance(DWF, PD, mean=prx_mean)
                prx_std = np.sqrt(prx_var)
                PEst.get_variance(DWF, PD, mean=PEst.mean_value)
                # print("Estimating proxy quartiles ...")
                prx_pc05 = get_proxy_percentile(5., cdfmat, PDist.proxyspan)
                prx_pc10 = get_proxy_percentile(10., cdfmat, PDist.proxyspan)
                prx_pc25 = get_proxy_percentile(25., cdfmat, PDist.proxyspan)
                prx_pc50 = get_proxy_percentile(50., cdfmat, PDist.proxyspan)
                prx_pc75 = get_proxy_percentile(75., cdfmat, PDist.proxyspan)
                prx_pc90 = get_proxy_percentile(90., cdfmat, PDist.proxyspan)
                prx_pc95 = get_proxy_percentile(95., cdfmat, PDist.proxyspan)

                
                baseRow = 5
                baseCol = 10
                colName = ['calage','mean','std','pc50','pc05','pc10','pc25','pc75','pc90','pc95']
                colVars = [calage, prx_mean, prx_std, prx_pc50, prx_pc05,prx_pc10, prx_pc25, prx_pc75,prx_pc90,prx_pc95]
                colLen = len(colName)
                for idx3 in range(colLen):
                    # break
                    colVar = colVars[idx3]
                    for idx4, val4 in enumerate(colVar):
                        if idx4==0:
                            if idx3==0:
                                outWSs[0].cell(baseRow-1,baseCol+iTarget*(colLen+2)+idx3).value=name2
                            outWSs[0].cell(baseRow+idx4,baseCol+iTarget*(colLen+2)+idx3).value=colName[idx3]
                        outWSs[0].cell(baseRow+1+idx4,baseCol+iTarget*(colLen+2)+idx3).value=val4
                # Current pointer
                if len(dataAll[iSheet][idx1][idx2])<5:
                    dataAll[iSheet][idx1][idx2]
                dataAll[iSheet][idx1][idx2][3]=copy.deepcopy(PDist)
                dataAll[iSheet][idx1][idx2][4]=copy.deepcopy(PEst)
#%%
        # wsheets[iSheet].cell(addresso[idx5][0]+1+idx4, addresso[idx5][1]).value = X1[idx5]    

    if bErrorCheck and bError:
        exit()

    with open('ProxyRecord.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump(dataAll, f)

else:
    # If the proxy record calculation was done already,
    print('Reading Preprocessed data (ProxyRecord.pkl)')
    for name in nameSheets:
        if name not in wb.sheetnames:
            outWSs.append(wb.create_sheet(name))
        else:
            outWSs.append(wb[name])

    # Read the dataAll from ProxyRecord.pkl
    with open('ProxyRecord.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        dataAll = pickle.load(f)

# wb.save(xlName)



# Obtain two trace element data from the Excel sheet
# Assumes the efficiecy is 1
K_e1, K_e2 = 1, 1 # Ni 1, Co 2
# Obtain the ln kinetic rates and standard deviations of trace elements
Kd_mn1, Kd_sd1 = wsheets[1].cell(6,5).value, wsheets[1].cell(6,6).value
Kd_mn2, Kd_sd2 = wsheets[1].cell(8,5).value, wsheets[1].cell(8,6).value
# 
# First trace element settings
TE1 = {'elem':wsheets[1].cell(6,4).value, 'mol_wt':wsheets[1].cell(6,7).value, 
       'Kp':wsheets[1].cell(6,8).value, 'temp_C':wsheets[1].cell(10,2).value,
       'aq_conc':np.float64(wsheets[1].cell(6,9).value), 
       'ca_conc':np.float64(wsheets[1].cell(11,2).value),
       'PDist':copy.deepcopy(dataAll[1][0][0][3])}
# Second trace element settings
TE2 = {'elem':wsheets[1].cell(8,4).value, 'mol_wt':wsheets[1].cell(8,7).value, 
       'Kp':wsheets[1].cell(8,8).value, 'temp_C':wsheets[1].cell(10,2).value,
       'aq_conc':np.float64(wsheets[1].cell(8,9).value), 
       'ca_conc':np.float64(wsheets[1].cell(11,2).value),
       'PDist':copy.deepcopy(dataAll[1][1][0][3])}

# Load the probability distribution of driprate (drips/minute) vs. time
# using the first trace element
# Span is the span of the drip rates
V_pdf_TE1, V_age_TE1, V_span_TE1 = driprates(Kd_mn1, Kd_sd1, K_e1,
                                            TE=TE1, calib=False, pbar_on=False)
# Load the probability distribution of driprate (drips/minute) vs. time
# using the second trace element
V_pdf_TE2, V_age_TE2, V_span_TE2 = driprates(Kd_mn2, Kd_sd2, K_e2,
                                            TE=TE2, calib=False, pbar_on=False)

# Check whether the age and the span of the drip rates are identical
assert np.all(V_age_TE2 == V_age_TE1), "ages do not match"
assert np.all(V_span_TE2 == V_span_TE1), "spans do not match"
V_age = V_age_TE1
V_span = V_span_TE1
# Construct span values for numerical integration
V_rsw = 0.5 * np.r_[
                    V_span[1] - V_span[0],
                    V_span[2:] - V_span[:-2],
                    V_span[-1] - V_span[-2]
                    ]

# take average of the two
# V_pdf = 0.5 * (V_pdf_TE2 + V_pdf_TE1)
# Construct a joint probability distribution of the drip rates of two trace elements
prod = V_pdf_TE1 * V_pdf_TE2
i = prod < 0.               # these are infinitesimally tiny numbers e-300
prod[i] = 0.
V_pdf = np.sqrt(prod)
C = (V_pdf.T * V_rsw).sum(axis=1)
# Normalize the distribution to make it pdf
V_pdf = V_pdf / C

# Evaluate the median drip rates from the joint pdf
# V_med = np.zeros(len(V_age))
# for i in range(V_pdf.shape[1]):
#     cdf = np.cumsum(V_rsw * V_pdf[:, i])
#     f_ = interp1d(cdf, V_span, kind="linear",
#                     bounds_error=False, fill_value=(0., 1.))
#     V_med[i] = f_(0.5)

# Evaluate the drip rates from the joint pdf
# baseRow and baseCol are starting row and column of the drip rate table
baseRow = 5
baseCol = 1
V_med = np.zeros(len(V_age))
V_pc05 = np.zeros(len(V_age))
V_pc10 = np.zeros(len(V_age))
V_pc25 = np.zeros(len(V_age))
V_pc75 = np.zeros(len(V_age))
V_pc90 = np.zeros(len(V_age))
V_pc95 = np.zeros(len(V_age))
for i in range(V_pdf.shape[1]):
    cdf = np.cumsum(V_rsw * V_pdf[:, i])
    f_ = interp1d(cdf, V_span, kind="linear",
                    bounds_error=False, fill_value=(0., 1.))
    # if idx4==0:
    #     if idx3==0:
    #         outWSs[0].cell(baseRow-1,baseCol+iTarget*(colLen+2)+idx3).value=name2
    #     outWSs[0].cell(baseRow+idx4,baseCol+iTarget*(colLen+2)+idx3).value=colName[idx3]
    V_med[i] = f_(0.5)
    V_pc05[i] = f_(0.05)
    V_pc10[i] = f_(0.10)
    V_pc25[i] = f_(0.25)
    V_pc75[i] = f_(0.75)
    V_pc90[i] = f_(0.90)
    V_pc95[i] = f_(0.95)
    # Set the values in the Excel sheet
    if i==0:
        for idx3, it1 in enumerate(['age','pc05','pc10','pc25','med','pc75','pc90','pc95']):
            if idx3==0: outWSs[1].cell(baseRow-1,baseCol+idx3).value="Probabilistically estimated range of drip rate vs. Time"
            outWSs[1].cell(baseRow,baseCol+idx3).value=it1
    outWSs[1].cell(baseRow+1+i,baseCol).value=V_age[i]
    outWSs[1].cell(baseRow+1+i,baseCol+1).value=V_pc05[i]
    outWSs[1].cell(baseRow+1+i,baseCol+2).value=V_pc10[i]
    outWSs[1].cell(baseRow+1+i,baseCol+3).value=V_pc25[i]
    outWSs[1].cell(baseRow+1+i,baseCol+4).value=V_med[i]
    outWSs[1].cell(baseRow+1+i,baseCol+5).value=V_pc75[i]
    outWSs[1].cell(baseRow+1+i,baseCol+6).value=V_pc90[i]
    outWSs[1].cell(baseRow+1+i,baseCol+7).value=V_pc95[i]


# Plot symbol and text size settings
majtiksz = MAJTIKSZ + 3
mintiksz = MINTIKSZ + 2
tiklabfs = TIKLABFS + 2
axlabfs = AXLABFS + 4
tiklabfs = TIKLABFS + 3

# Plot drip rate vs time
# First plot configuration
print('Creating drip rate vs. time plot')
plotConf = plotConfs[0]
fig = pl.figure(figsize=[plotConf['Width'], plotConf['Height']])
ax1 = fig.add_axes([0.10, 0.125, 0.85, 0.80])
V_ = V_pdf / np.max(V_pdf, axis=0)
im1 = pl.pcolormesh(V_age, V_span, V_, cmap=pl.cm.Blues, rasterized=True)

# im1 = pl.pcolormesh(V_age, V_span, V_, cmap=pl.cm.Blues, rasterized=True)
cb1 = pl.colorbar(im1)
# clr = "Tomato"
# pd.DataFrame({'age':V_age,'V_med':V_med,'V_pc25':V_pc25,'V_pc75':V_pc75}).to_clipboard(index=False,header=False)
ax1.plot(V_age, V_med, plotConf['Plot1Type'], c=plotConf['Plot1Clr'], label=plotConf['Plot1Label'])
ax1.plot(V_age, V_pc25, plotConf['Plot2Type'], c=plotConf['Plot2Clr'], label=plotConf['Plot2Label'])
ax1.plot(V_age, V_pc75, plotConf['Plot3Type'], c=plotConf['Plot3Clr'], label=plotConf['Plot3Label'])
ax1.set_xlabel(plotConf['Xlabel'], fontsize=axlabfs)
ax1.set_ylabel(plotConf['Ylabel'], fontsize=axlabfs)
ax1.set_xlim(plotConf['Xmin'], plotConf['Xmax'])
ax1.xaxis.set_major_locator(ticker.MultipleLocator(plotConf['XtickMajor'])) # set BIG ticks
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(plotConf['XtickMinor'])) # set BIG ticks
ax1.yaxis.set_major_locator(ticker.MultipleLocator(plotConf['YtickMajor'])) # set BIG ticks
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(plotConf['YtickMinor'])) # set BIG ticks
# ax1.set_xlim(ax1.get_xlim()[::-1])
ax1.set_ylim(plotConf['Ymin'], plotConf['Ymax'])
ax1.tick_params(labelsize=tiklabfs, size=majtiksz)
ax1.tick_params(which="minor", size=mintiksz)
# ax1.set_xticks(np.arange(10000, -101, -1000))
# ax1.set_xticks(np.arange(plotConf['Xmin'], plotConf['Xmax'], plotConf['Xtick']), minor=True)
# ax1.set_yticks(np.arange(0., 35.1, 5))
# ax1.set_yticks(np.arange(plotConf['Ymin'], plotConf['Ymax'], plotConf['Ytick']), minor=True)
ax1.grid(which="both",alpha=0.25)
ax1.legend(loc="lower left")
cb1.set_label(plotConf['CMapLabel'], fontsize=axlabfs)
cb1.ax.tick_params(labelsize=tiklabfs, size=majtiksz)
cb1.ax.tick_params(which="minor", size=mintiksz)


outPlot = outWSs[3]
for _ in enumerate(outPlot._images):
    del(outPlot._images[0])
outPlot.cell(1,1).value = 'Drip rate (per min) vs. Time (BP)'
pltMain = io.BytesIO()
pl.savefig(pltMain, format='png', dpi=200, bbox_inches='tight')
pltMain.seek(0)
im = Image.open(pltMain)
img = openpyxl.drawing.image.Image(im)
img.anchor = 'A3'
outPlot.add_image(img)
# pltMain.close()


# Plot drip rate and isotope
if plotConfs[1]['Plot2Isotope'] in [1,2]:
    print('Creating drip rate and isotope plot')
    # Calculate Isotope median value
    PDist1 = copy.deepcopy(dataAll[2][plotConfs[1]['Plot2Isotope']-1][0][3])
    cdfmat = get_cdfmat(PDist1.pdfmat, PDist1.proxyspan)
    prx_pc50 = get_proxy_percentile(50., cdfmat, PDist1.proxyspan)
    for i in range(V_pdf.shape[1]):
        if i==0:
            for idx3, it1 in enumerate(['age','median']):
                if idx3==0: outWSs[2].cell(baseRow-1,baseCol+idx3).value="Probabilistically estimated median of isotope data vs. Time"
                outWSs[2].cell(baseRow,baseCol+idx3).value=it1
        outWSs[2].cell(baseRow+1+i,baseCol).value=V_age[i]
        outWSs[2].cell(baseRow+1+i,baseCol+1).value=prx_pc50[i]

    plotConf = plotConfs[1]
    fig = pl.figure(figsize=[plotConf['Width'], plotConf['Height']])
    ax1 = fig.add_axes([0.10, 0.125, 0.85, 0.80])

    # pd.DataFrame({'age':V_age,'V_med':V_med,'V_pc25':V_pc25,'V_pc75':V_pc75}).to_clipboard(index=False,header=False)
    pl1 = ax1.plot(V_age, V_med, plotConf['Plot1Type'], c=plotConf['Plot1Clr'], label=plotConf['Plot1Label'])
    ax1.set_xlabel(plotConf['Xlabel'], fontsize=axlabfs)
    ax1.set_ylabel(plotConf['Y1label'], fontsize=axlabfs)
    ax1.set_xlim(plotConf['Xmin'], plotConf['Xmax'])
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(plotConf['XtickMajor'])) # set BIG ticks
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(plotConf['XtickMinor'])) # set BIG ticks
    # ax1.set_xlim(ax1.get_xlim()[::-1])
    ax1.set_ylim(plotConf['Y1min'], plotConf['Y1max'])
    ax1.tick_params(labelsize=tiklabfs, size=majtiksz)
    ax1.tick_params(which="minor", size=mintiksz)
    # ax1.set_xticks(np.arange(10000, -101, -1000))
    # ax1.set_xticks(np.arange(-plotConfs[0]['Xmin'], -plotConfs[0]['Xmax'],-plotConfs[0]['Xtick']), minor=True)
    # ax1.set_yticks(np.arange(-20, 0, 1))
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(plotConf['Y1tickMajor'])) # set BIG ticks
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(plotConf['Y1tickMinor'])) # set BIG ticks
    # ax1.set_yticks(np.arange(plotConfs[0]['Ymin'], plotConfs[0]['Ymax'], plotConfs[0]['Ytick']), minor=True)
    ax1.grid(which="both",alpha=0.25)
    ax1.legend(loc="lower left")
    # secax = ax1.secondary_yaxis('right')
    # ax2 = fig.add_axes([0.10, 0.125, 0.85, 0.80])
    ax2 = ax1.twinx()
    pl2 = ax2.plot(V_age, prx_pc50, plotConf['Plot2Type'], c=plotConf['Plot2Clr'], label=plotConf['Plot2Label'])
    ax2.set_ylim(plotConf['Y2min'], plotConf['Y2max'])
    # ax2.set_ylim(-6, -12)
    ax2.set_ylabel(plotConf['Y2label'], fontsize=axlabfs)
    lns = pl1+pl2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc="lower left")
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(plotConf['Y2tickMajor'])) # set BIG ticks
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(plotConf['Y2tickMinor'])) # set BIG ticks

    outPlotIso = outWSs[4]
    for _ in enumerate(outPlotIso._images):
        del(outPlotIso._images[0])
    outPlotIso.cell(1,1).value = 'Drip rate (per min) vs. Isotope signature'
    pltMain = io.BytesIO()
    pl.savefig(pltMain, format='png', dpi=200, bbox_inches='tight')
    pltMain.seek(0)
    im = Image.open(pltMain)
    img = openpyxl.drawing.image.Image(im)
    img.anchor = 'A3'
    outPlotIso.add_image(img)

# Save the Excel workbook
wb.save(xlName)
