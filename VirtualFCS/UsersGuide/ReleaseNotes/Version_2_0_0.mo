class Version_2_0_0
  extends VirtualFCS.UsersGuide.ReleaseNotes;
  
  
  annotation(Documentation(info = "<html></p>

<p><head> <meta http-equiv=Content-Type content=&quot;text/html; charset=windows-1252&quot;> <meta name=Generator content=&quot;Microsoft Word 15 (filtered)&quot;> <style> <!--  /* Font Definitions */  @font-face  {font-family:Wingdings;  panose-1:5 0 0 0 0 0 0 0 0 0;} @font-face  {font-family:&quot;Cambria Math&quot;;  panose-1:2 4 5 3 5 4 6 3 2 4;} @font-face  {font-family:Calibri;  panose-1:2 15 5 2 2 2 4 3 2 4;}  /* Style Definitions */  p.MsoNormal, li.MsoNormal, div.MsoNormal  {margin-top:0cm;  margin-right:0cm;  margin-bottom:8.0pt;  margin-left:0cm;  line-height:107%;  font-size:11.0pt;  font-family:&quot;Calibri&quot;,sans-serif;} p.MsoHeader, li.MsoHeader, div.MsoHeader  {mso-style-link:&quot;Header Char&quot;;  margin:0cm;  font-size:11.0pt;  font-family:&quot;Calibri&quot;,sans-serif;} a:link, span.MsoHyperlink  {color:#0563C1;  text-decoration:underline;} p.MsoListParagraph, li.MsoListParagraph, div.MsoListParagraph  {margin-top:0cm;  margin-right:0cm;  margin-bottom:8.0pt;  margin-left:36.0pt;  line-height:107%;  font-size:11.0pt;  font-family:&quot;Calibri&quot;,sans-serif;} p.MsoListParagraphCxSpFirst, li.MsoListParagraphCxSpFirst, div.MsoListParagraphCxSpFirst  {margin-top:0cm;  margin-right:0cm;  margin-bottom:0cm;  margin-left:36.0pt;  line-height:107%;  font-size:11.0pt;  font-family:&quot;Calibri&quot;,sans-serif;} p.MsoListParagraphCxSpMiddle, li.MsoListParagraphCxSpMiddle, div.MsoListParagraphCxSpMiddle  {margin-top:0cm;  margin-right:0cm;  margin-bottom:0cm;  margin-left:36.0pt;  line-height:107%;  font-size:11.0pt;  font-family:&quot;Calibri&quot;,sans-serif;} p.MsoListParagraphCxSpLast, li.MsoListParagraphCxSpLast, div.MsoListParagraphCxSpLast  {margin-top:0cm;  margin-right:0cm;  margin-bottom:8.0pt;  margin-left:36.0pt;  line-height:107%;  font-size:11.0pt;  font-family:&quot;Calibri&quot;,sans-serif;} span.HeaderChar  {mso-style-name:&quot;Header Char&quot;;  mso-style-link:Header;} .MsoChpDefault  {font-family:&quot;Calibri&quot;,sans-serif;} .MsoPapDefault  {margin-bottom:8.0pt;  line-height:107%;}  /* Page Definitions */  @page WordSection1  {size:595.3pt 841.9pt;  margin:72.0pt 72.0pt 72.0pt 72.0pt;} div.WordSection1  {page:WordSection1;}  /* List Definitions */  ol  {margin-bottom:0cm;} ul  {margin-bottom:0cm;} --> </style></p>

<p></head></p>

<p><body lang=EN-GB link=&quot;#0563C1&quot; vlink=&quot;#954F72&quot; style='word-wrap:break-word'></p>

<p><div class=WordSection1></p>

<p><p class=MsoNormal><span style='display:none'>Top of Form</span></p></p>

<p><ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Introduction</li>  <ul style='margin-top:0cm' type=disc>  <li class=MsoNormal>Brief overview of the major upgrade release</li>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Major changes include, </li>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Library is upgraded to be compatible with OpenModelica v.1.20.  </li>  <li class=MsoNormal>A new package, ComponentTesting, is introduced to the  library. </li>  <li class=MsoNormal>Computational time reduced. </li>  </ol>  <li class=MsoNormal>Minor changes include,</li>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Changed bitmaps to vector graphics.</li>  <li class=MsoNormal>Restructuring of library following the Modelica  structure.</li>  <li class=MsoNormal>Rolling force was re-introduced in the Vehicle model.</li>  <li class=MsoNormal>Removed i<sup>2</sup>R terms from heat sources in the thermal  part for Fuel Cell and Battery models.</li>  </ol>  </ol>  <li class=MsoNormal>List of key features and improvements included in the  upgrade,</li>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Performance indicating parameters such power and various  efficiencies for vehicle, fuel cell stack, battery, and BOP components  are introduced.</li>  <li class=MsoNormal>Introduced dynamic control of BOP components such as Air  subsystem, Hydrogen subsystem and Cooling subsystem. </li>  <li class=MsoNormal>Faster simulations of fluids library from previous  version. </li>  </ol>  </ul>  <li class=MsoNormal>New Features</li>  <ul style='margin-top:0cm' type=disc>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Library upgrade: Virtual-FCS is upgraded to be  compatible with OpenModelica v.1.20 and the Modelica Standard Library  4.0.0 (previously only up to OM v.1.14 and MSL 3.2.3). The three balance  of plant subsystems, air, hydrogen, and cooling, all including fluid  components such as pipes, valves, tanks, sinks, and pumps have been  redeclared and partially remodelled to overcome compatibility issues  with the latest Modelica fluid library. The medium packages are initialised  and constrained more clearly which also enables designers to easily  change initial conditions of models.</li>  <li class=MsoNormal>New class: A new package, ComponentTesting, is  introduced to enable designers to debug larger models by testing simple models  for single components. </li>  <li class=MsoNormal>Faster simulations are available due to faster default  solver and library upgrades. </li>  <li class=MsoNormal>Changed bitmaps to vector graphic to avoid long and  unreadable model files and any compilation bugs related to bitmap. </li>  <li class=MsoNormal>The library is restructured complying to the Modelica Standard  library structure providing familiarity to users. </li>  <li class=MsoNormal>Calculations for power and efficiency for vehicle, fuel  cell stack, battery, and BOP components is introduced to help users  evaluate key performance indications. </li>  <li class=MsoNormal>Control of BOP components is now developed to  dynamically respond to varying fuel cell load. This includes adding  separate control blocks for different subsystem, which also enables  users to change only parts of the system.</li>  </ol>  </ul>  <li class=MsoNormal>Bug Fixes</li>  <ul style='margin-top:0cm' type=disc>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Modelica standard library&rsquo;s fluid system model is required  in the VirtualFCS library to use fluid components. In the previous  version this was introduced at every level in Virtual-FCS which was  redundant. This issue was resolved by only declaring the system at the  top level by using the inner keyword and referring to the same system in  the lower level by using the outer keyword. This improves stability of  models and avoids confusion for users. </li>  <li class=MsoNormal>Fuel cell and Battery were never reaching the desired/expected  temperature. This was fixed by altering the thermal parameters and change  the cooling control to be dynamic. </li>  <li class=MsoNormal>Battery size in few models had incorrect voltage. Fixed  this by changing the voltage for BOP components to be 24 V. </li>  <li class=MsoNormal>Selecting/deselecting to use regenerative breaking did  not function. This was fixed. </li>  <li class=MsoNormal>Converted few classes to models to resolve OM-warnings. </li>  </ol>  </ul>  <li class=MsoNormal>Deprecated Features</li>  <ul style='margin-top:0cm' type=disc>  <li class=MsoNormal>Library will now only support OpenModelica version 1.20 </li>  <li class=MsoNormal>Library will now only support Modelica Device Drivers  version 2.1.1</li>  </ul>  <li class=MsoNormal>Breaking Changes</li>  <ul style='margin-top:0cm' type=disc>  <li class=MsoNormal>With this upgrade one cannot use OM older versions</li>  </ul>  <li class=MsoNormal>Upgrade Instructions</li>  <ul style='margin-top:0cm' type=disc>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Install OpenModelica v 1.20 from <a  href=&quot;https://build.openmodelica.org/omc/builds/windows/releases/1.20/0/&quot;>https://build.openmodelica.org/omc/builds/windows/releases/1.20/0/</a>  </li>  <li class=MsoNormal>Open the User interface<span style='font-family:Wingdings'>&agrave;</span>  Install Library<span style='font-family:Wingdings'>&agrave;</span> Modelica Device  Driver v 2.1.1</li>  <li class=MsoNormal>Download/Get latest code from <a  href=&quot;https://github.com/Virtual-FCS/VirtualFCS&quot;>https://github.com/Virtual-FCS/VirtualFCS</a></li>  <li class=MsoNormal>Load the VirtualFCS library.</li>  </ol>  </ul>  <li class=MsoNormal>Support and Feedback</li>  <ul style='margin-top:0cm' type=disc>  <li class=MsoNormal>Users are encouraged to report their issues on <a  href=&quot;https://github.com/Virtual-FCS/VirtualFCS/issues&quot;>https://github.com/Virtual-FCS/VirtualFCS/issues</a>  if they seek support or feedback. </li>  </ul>  <li class=MsoNormal>Contributors</li>  <ul style='margin-top:0cm' type=disc>  <li class=MsoNormal>List of contributors made the development of this  release.</li>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Benjamin Synnev&aring;g</li>  <li class=MsoNormal>Mike Gerhardt</li>  </ol>  <li class=MsoNormal>List of contributors who helped on the development of  this release.</li>  <ol style='margin-top:0cm' start=1 type=1>  <li class=MsoNormal>Knut Vidar Skjersli</li>  <li class=MsoNormal>Yash Raka</li>  </ol>  </ul>  <li class=MsoNormal>Acknowledgements</li>  <ul style='margin-top:0cm' type=disc>  <li class=MsoNormal>Dietmar Winkler helped on the development of this release.</li>  </ul>  <li class=MsoNormal>Release date</li> </ol></p> <ul style='margin-top:0cm' type=disc>  <li class=MsoNormal>On date 19 January 2023 release was made available to the public.</li>  </ul>

<p></div></p>

<p></body></p>

<p></html>"));
end Version_2_0_0;
