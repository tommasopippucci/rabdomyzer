#!/usr/bin/env python
import sys
import argparse
import pandas as pd

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Create an excel worksheet with rabdomyzer results')
  parser.add_argument('--het', help="het file created with rabdomyzer", type=argparse.FileType('r'))
  parser.add_argument('--hom', help="hom file created with rabdomyzer", type=argparse.FileType('r'))
  parser.add_argument('--comhet', help="comhet file created with rabdomyzer", type=argparse.FileType('r'))
  parser.add_argument('--output', help="output file", type=str)
  if len(sys.argv)==1:
     print
     print "Usage: python excel_rabdomyzer.py --help"
     print
     raise StandardError("No arguments provided")
  arguments = parser.parse_args()

file_het= pd.read_csv(arguments.het, sep='\t',index_col = False, encoding='utf-8')
file_hom=pd.read_csv(arguments.hom, sep='\t', index_col = False, encoding='utf-8')
#file_het= pd.read_csv(arguments.het, sep='\t',index_col = False)
#file_hom=pd.read_csv(arguments.hom, sep='\t', index_col = False)
writer= pd.ExcelWriter(arguments.output, engine='xlsxwriter')
file_het.to_excel(writer, sheet_name='het_calls')
file_hom.to_excel(writer, sheet_name='hom_calls')
if len(sys.argv)>7:
   file_comhet=pd.read_csv(arguments.comhet, sep='\t', index_col = False, encoding='utf-8')
   #file_comhet=pd.read_csv(arguments.comhet, sep='\t', index_col = False)
   file_comhet.to_excel(writer, sheet_name='comhet_calls')
workbook  = writer.book
worksheet_het = writer.sheets['het_calls']
worksheet_hom = writer.sheets['hom_calls']

###########################Colors to use################################

format_red = workbook.add_format({'bg_color' : '#FF3333'})
format_green = workbook.add_format({'bg_color' : '#66FF66'})
format_orange = workbook.add_format({'bg_color' : '#FFA500'})

###########add conditional formatting to het file########################

worksheet_het.conditional_format('AK2:AK100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'tolerated',
                                       'format':   format_green})

worksheet_het.conditional_format('AK2:AK100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'deleterious',
                                       'format':   format_red})


worksheet_het.conditional_format('AL2:AL100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'damaging',
                                       'format':   format_red})

worksheet_het.conditional_format('AL2:AL100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'benign',
                                       'format':   format_green})

worksheet_het.conditional_format('AN2:AN100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'D',
                                       'format':   format_red})
                                       
worksheet_het.conditional_format('AN2:AN100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'T',
                                       'format':   format_green})
                                      
worksheet_het.conditional_format('AM2:AM100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '15',
                                       'format':   format_red})
                                       
worksheet_het.conditional_format('AO2:AO100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '3',
                                       'format':   format_red})                                                                              

worksheet_het.conditional_format('AP2:AP100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '0.5',
                                       'format':   format_red}) 
                                       
worksheet_het.conditional_format('AQ2:AQ100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '0.5',
                                       'format':   format_red})

############add conditional formatting to Hom_file##########################

worksheet_hom.conditional_format('AK2:AK100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'tolerated',
                                       'format':   format_green})

worksheet_hom.conditional_format('AK2:AK100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'deleterious',
                                       'format':   format_red})


worksheet_hom.conditional_format('AL2:AL100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'damaging',
                                       'format':   format_red})

worksheet_hom.conditional_format('AL2:AL100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'benign',
                                       'format':   format_green})

worksheet_hom.conditional_format('AN2:AN100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'D',
                                       'format':   format_red})
                                       
worksheet_hom.conditional_format('AN2:AN100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'T',
                                       'format':   format_green})
                                      
worksheet_hom.conditional_format('AM2:AM100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '15',
                                       'format':   format_red})
                                       
worksheet_hom.conditional_format('AO2:AO100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '3',
                                       'format':   format_red}) 

worksheet_hom.conditional_format('AP2:AP100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '0.5',
                                       'format':   format_red}) 
                                       
worksheet_hom.conditional_format('AQ2:AQ100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '0.5',
                                       'format':   format_red})
                                        
######################format comhet file if it exists####################

if len(sys.argv)>7:
    worksheet_comhet = writer.sheets['comhet_calls']
    
    worksheet_comhet.conditional_format('AK2:AK100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'tolerated',
                                       'format':   format_green})

    worksheet_comhet.conditional_format('AK2:AK100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'deleterious',
                                       'format':   format_red})


    worksheet_comhet.conditional_format('AL2:AL100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'damaging',
                                       'format':   format_red})

    worksheet_comhet.conditional_format('AL2:AL100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'benign',
                                       'format':   format_green})

    worksheet_comhet.conditional_format('AN2:AN100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'D',
                                       'format':   format_red})
                                       
    worksheet_comhet.conditional_format('AN2:AN100000', {'type':     'text',
                                       'criteria': 'containing',
                                       'value':    'T',
                                       'format':   format_green})
                                      
    worksheet_comhet.conditional_format('AM2:AM100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '15',
                                       'format':   format_red})
                                       
    worksheet_comhet.conditional_format('AO2:AO100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '3',
                                       'format':   format_red})
                                       
    worksheet_comhet.conditional_format('AP2:AP100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '0.5',
                                       'format':   format_red}) 
                                       
    worksheet_comhet.conditional_format('AQ2:AQ100000', {'type':     'cell',
                                       'criteria': 'greater than or equal to',
                                       'value':    '0.5',
                                       'format':   format_red}) 

writer.save()
