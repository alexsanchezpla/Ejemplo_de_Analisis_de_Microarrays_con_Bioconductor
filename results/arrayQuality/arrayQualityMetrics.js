// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, true, false, false, false, false, false, false, true, true, true, false, false ];
var arrayMetadata    = [ [ "1", "GSM26878.CEL", "GSM26878", "PF14 EnPnT2N1G2", "PF14", "A", "3", "18", "07/17/03 07:05:19" ], [ "2", "GSM26883.CEL", "GSM26883", "PF19 EnPuT4N0Gu", "PF19", "A", "2", "18", "06/11/03 13:41:28" ], [ "3", "GSM26887.CEL", "GSM26887", "PF23 EnPnT2N0G2", "PF23", "A", "3", "18", "07/22/03 05:45:18" ], [ "4", "GSM26903.CEL", "GSM26903", "PF39 EnPuT4N0Gu", "PF39", "A", "3", "17", "07/17/03 06:28:12" ], [ "5", "GSM26910.CEL", "GSM26910", "PF46 EnPnT4N1G3", "PF46", "A", "3", "19", "08/07/03 10:16:11" ], [ "6", "GSM26888.CEL", "GSM26888", "PF24 EnPnTiN0G3", "PF24", "B", "2", "NA", "07/17/03 05:22:10" ], [ "7", "GSM26889.CEL", "GSM26889", "PF25 EnPnT3N2G2", "PF25", "B", "1", "21", "07/18/03 04:51:27" ], [ "8", "GSM26892.CEL", "GSM26892", "PF28 EnPnT2N1G3", "PF28", "B", "2", "22", "07/17/03 04:53:32" ], [ "9", "GSM26898.CEL", "GSM26898", "PF34 EnPnT3N1G3", "PF34", "B", "1", "17", "07/18/03 06:36:27" ], [ "10", "GSM26906.CEL", "GSM26906", "PF42 EnPnT2N2G3", "PF42", "B", "na", "16", "07/18/03 04:39:49" ], [ "11", "GSM26879.CEL", "GSM26879", "PF15 EpPnTiN1G3", "PF15", "L", "2", "16", "07/22/03 04:19:40" ], [ "12", "GSM26896.CEL", "GSM26896", "PF32 EpPnT3N1G2", "PF32", "L", "2", "23", "07/22/03 04:59:24" ], [ "13", "GSM26897.CEL", "GSM26897", "PF33 EpPnTiN0G2", "PF33", "L", "1", "18", "07/17/03 04:59:21" ], [ "14", "GSM26907.CEL", "GSM26907", "PF43 EpPpT2N1G2", "PF43", "L", "1", "20", "08/07/03 10:50:30" ], [ "15", "GSM26911.CEL", "GSM26911", "PF47 EpPpT3N1G3", "PF47", "L", "3", "17", "08/07/03 10:44:53" ] ];
var svgObjectNames   = [ "pca", "dens", "dig" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
	var success = false;
	i = 0; 
	/* Some of this looping could already be cached in reportInit() */
	while( (!success) & (i < ssrules.length) ) {
	    selector = ssrules[i].selectorText;  // The selector 
            if (!selector) 
		continue; // Skip @import and other nonstyle rules
            if (selector == (".aqm" + reportObjId)) {
		success = true; 
		ssrules[i].style.cssText = cssText[0+status];
	    } else {
		i++;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
