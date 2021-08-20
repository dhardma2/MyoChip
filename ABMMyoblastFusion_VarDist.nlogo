extensions[matrix]
turtles-own [prolif? prolif-tick res-tick res-tot myotube? neighbours tubeset? tubeno d dist xa ya xnew ynew xvect yvect fusion-age MT-fuseable?]
globals [mt chance filenumber fusecount tickno tc tube-reset tn  x y linreg c m i  xcent ycent xb yb x1 y1 maximum j global-linkdist-list coeff-var coeff-var-list global-link-count data-collect? linkdistlist largelink
xminus? xplus? yminus? yplus?]
patches-own []
breed[myocytes myocyte]
breed[myotubes myotube]
myotubes-own[x_2 y_2]
links-own[largest? force-check?]
to Rec
  ;for recording images of frames
  set filenumber 1
  repeat 100
  [ go
    export-view (word "image" filenumber ".png")
    set filenumber filenumber + 1 ]
end

to setup
  clear-all
  set data-collect? false
  set global-linkdist-list []
  set coeff-var-list []
  set coeff-var 0
  set global-link-count []
  set linkdistlist []
  ;Initialise variables and arrays
  set tickno 0
  set fusecount 0
  set x []
  set y []
  set i 1
  set linreg []
  set tube-reset 0
  ;seed initial combined nuclei from count given by 'myoblast-seed' slider
  create-myocytes myoblast-seed
  [ ;initialise nuclei variables
    setxy random-xcor random-ycor
    myocreate]

  set-default-shape turtles "circle"
  reset-ticks
  ask turtles
    [

      ;set proportion of nuclei which belong to myotubes given by 'myotube-seed' slider (proportion is chosen stochastically)
      if (random 100 < Myotube-seed)
      [
        set breed myotubes
        set res-tot 0
        set x_2 xcor
        set y_2 ycor
    ]
  ]
end

to go
  clear-links
  ;option to show the residence time ticks for all nuclei
  ;showresticks
  if ticks > 7200
  [
    set data-collect? true
  ]

  ask myocytes
  [
    ;if the nucleus is not currently proliferating
    ifelse prolif-tick = 0
      [
        ;change angle of motion with degree of turn chosen from a normal distribution with mean given by the 'angle' slider and standard deviation by the 'anglesd' slider
        rt random-normal angle anglesd
        ;move forward with distance chosen from a normal distribution with mean given by the 'speed' slider and standard deviation by the 'speedsd' slider
        forward random-normal speed speedsd
      ]
    ;if the nucleus is currently proliferating, increase the proliferation time by 1 minute.
      [set prolif-tick prolif-tick + 1]
    ;act of proliferation takes 20 minutes
if prolif-tick > 20
    [ set prolif-tick 0
     set res-tot 0]
  ]
  ask turtles
  [
  ;myocyte nuclei within the radius of influence of a myotube nucleus get an uptick in their residence time
  restime

    ;setting a maximum total residence time of 1 per tick removes the effects of interactions of clusters of nuclei in ramping up residence time
    ifelse res-tick > 0
    [
        set res-tot res-tot + 1
        if res-tot > resmax + 10
        [set res-tot resmax + 10]
    ]
     [set res-tot 0]
  ]
  ;do not begin proliferation until the residence time threshold is reached to allow comparison of models at T0 = Tresmax
  if tickno >= resmax
   [proliferate]

ask myocytes
  [
    let RM1 random-normal resmax resmaxsd
    if RM1 <= 0
    [set RM1 1]
  ;once the residence time surpasses the total, myocytes fuse.
  if res-tot >= RM1
    [myogen]
  ]
if tube-reset > 0
  [
    ask myotubes
    [
      if tubeno = tube-reset
      [set res-tot 0]
    ]
  ]
ask myotubes
  [
    ;If myotube nuclei are younger than threshold age for fusion, set myotube nuclei into discrete myotubes with an index number
    if  MT-fuseable? = true
    [
      tube-define
    ]
    ;add one minute to fusion age
    set fusion-age fusion-age + 1
    let MMTA Max-MT-age
    if fusion-age > (MMTA * 1440 + resmax)
    [set MT-fuseable? false]

  ]
  set maximum max [tubeno] of myotubes
  if ticks > 1
 [
    if Nuclei_force
    ;apply lateral force to preserve tube shape
    ;Apply nuclei repulsion force to prevent overlap and apply actin chain elongation
    [nuclei-force]
  ]
  ask myotubes
  [myomove]

  if data-collect? = true
  [
    make-links
    data-collect-global
  ]
  set data-collect? false
  tick
  set tickno tickno + 1
end

to myocreate
  set chance 1 / chance-division
  set prolif? false
  set prolif-tick 0
  set res-tick 0
  set res-tot 0
  set tubeset? false
  set myotube? false
  set fusion-age 0
  set MT-fuseable? true
  set xvect 0
  set yvect 0
  ;random distribution of cells
  if GUI
  [
    set color red
    set size  radius / 2
  ]
  set dist 0

end

to proliferate
   ;myocyte nuclei have a chance of proliferating at each timestep given by the 'chance-division' slider .
    ask myocytes
     [
      if prolif-tick = 0
      [
        if ((random chance) = 0)
        [
        set prolif? true
        set prolif-tick 1
        ;create a new myocyte and move it away at 180 deg to the original
        hatch 1 [lt 180]
     ]
    ]
   ]

end

to restime
  ;if there are unproliferating nuclei within the RoI then increase the neighbour's reidence time ticker by 1 (neighbour will do the same for you!)
  set neighbours other turtles in-radius (radius)
  ifelse any? neighbours
  [
   ask neighbours
    [
    if prolif-tick = 0 and res-tick = 0

    [set res-tick res-tick + 1]
    ]
  ]
   [set res-tick 0]
end

to myogen
  ;ensure that fusing myocyte has a neighbouring nuclei which is also ready to fuse
  set neighbours other turtles in-radius (radius)
  set xminus? false
  set xplus? false
  set yminus? false
  set yplus? false
  if any? neighbours
  [
   ask neighbours
    [
      let RM2 random-normal resmax resmaxsd
      if RM2 <= 0
      [set RM2 1]
    if res-tot >= RM2
      [
        ;neighbour becomes (or is already) a myotube. Residence time ticker resets.
        if breed = myocytes
        [
          set breed myotubes
          set x_2 xcor
          set y_2 ycor
        ]
       if breed = myotubes
        [

          if x_2 < min-pxcor
          [set xminus? true]
          if x_2 > max-pxcor
          [set xplus? true]
          if y_2 < min-pycor
          [set yminus? true]
          if y_2 > max-pycor
          [set yplus? true]

        ]


        set mt 1
        set res-tot 0
        ;option to allow only one nuclei to fuse at a time. All myotube nuclei within a cell
        ;have their residence time ticks reset to 0 and so it will take at least until t=resmax
        ;until another cell can fuse.
        if single-fusion
        [one-at-a-time]
      ]
    ]
    if mt = 1
    [
    ;change to myotube nuclei, add 1 to fuse count and reset residence time ticker.
      set breed myotubes
      set x_2 xcor
      set y_2 ycor
      if xminus? = true
      [ set x_2 x_2 - max-pxcor]
      if xplus? = true
      [set x_2 x_2 + max-pxcor]
      if yminus? = true
      [set y_2 y_2 - max-pycor]
      if yplus? = true
      [set y_2 y_2 + max-pycor]
      set fusecount fusecount + 1
      set mt 0
      set res-tot 0
    ]
  ]

end
to tube-define
  let myx_2 x_2
  let myy_2 y_2
  set neighbours other myotubes in-radius (radius)

  set tc 0
  set tn 0
  ;if the nuclei belongs to a numbered myotube, set any neighbouring nuclei to the same myotube.
  ifelse tubeset? = true
  [
    set tn tubeno
    if any? neighbours
    [
      ask neighbours
      [
        if (abs (myx_2 - x_2) < 100) and (abs (myy_2 - y_2) < 100)
      [
        set tubeset? true
        set tubeno tn
      ]
      ]
    ]
  ]
  ;if the nuclei does not belong to a numbered myotube, look for neighbours.
  [
    ifelse any? neighbours
      [
        ;if nuclei has neighbours belonging to a numbered myotube then add current nuclei to this myotube.
        ask neighbours
        [
          if (abs (myx_2 - x_2) < 100) and (abs (myy_2 - y_2) < 100)
          [
          ifelse tubeset? = true
          [
            set tn tubeno
          ]
          ;if nuclei has neighbours not belonging to a numbered myotube then define a new myotube and set the neighbours in it.
          [
            set tubeset? true
            ;set the myotube index number to the myocyte unique turtle number ('who')
            set tn who + 1
            set tubeno tn
          ]
          ]
        ]
    ]
    ;if nuclei has no neighbours (only possible at t=1) then then define a new myotube.
    [
      set tubeset? true
      ;set the myotube index number to the myocyte unique turtle number ('who')
      set tn who + 1
    ]
    ;Place current nuclei in newly defined myotube.
    ]  set tubeno tn

  set color red
end

to nuclei-force
  while [i <= maximum ]
  [
    ;create arrays of x and y coordinates of nuclei in each myotube
    set x1 [x_2] of myotubes with [tubeno = i]
    set y1 [y_2] of myotubes with [tubeno = i]
    ;for myotubes with > 2 nuclei...
    if length x1 > 2
    [
      ;find centroid of nuclei
      set xcent mean x1
      set ycent mean y1

      ask myotubes with [tubeno = i]
      [
        ;set axis at centroid **can this be done without looping?**
        set x lput (x_2 - xcent) x
        set y lput (y_2 - ycent) y
      ]
      ;apply least-squares linear regression to find gradient (m) of nuclei
      set linreg item 0 (matrix:regress matrix:from-column-list (list y x))
      set m item 1 linreg
      ;set points normal to centreline and move nuclei towards centreline
      ptmove
    ]
    set x[]
    set y[]
    set i i + 1
  ]
  set i 1
  ;calculate inter-nuclei repulsion force
  while [i <= maximum ]
  [
    ;for tubes with more than 2 nuclei
    if count myotubes with [tubeno = i] > 2
    [
    ask myotubes with [tubeno = i]
      [
        let whoself who
        let xi x_2
        let yi y_2
        let xsum 0
        let ysum 0
        let b 0
        let a 0
        ;cycle through other nuclei to calculate summative force
        ask myotubes with [tubeno = i]
        [
          if who != whoself
          [
            set a (xi - x_2)
            set b (yi - y_2)
            ;strategy for if repulsive force diminishes with distance from nuclei.
            let NF random-normal nuc-forceB nuc-forcesd
            ifelse DistDep = true
            [
              let force 0
              let csq (a ^ 2) + (b ^ 2)

              if csq != 0
              [
                ;f(nuc)= (d/nuc-forceB)^-1

                set force NF / csq
                ;apply bounds on force
                if abs force < 1.0E-5
                [
                  ;print "error:force out of bounds (low)"
                  ;print force
                  ;print i
                  set force 0
                ]
                if abs force > 1.5
                [
                  ;print "error:force out of bounds (high)"
                  ;print force
                  ;print i
                  ;set force 1.5
                ]
                set csq 0
              ]
              ;sum forces in x and y directions
              set xsum xsum + force * a
              set ysum ysum + force * b

            ]
            [
              if abs a > 0
              [set xsum xsum + (NF * (a / abs a))]
              if abs b > 0
              [set ysum ysum + (NF * (b / abs b))]

            ]
          ]
        ]
        ;apply bounds on sum of forces
        if abs xsum > 2
                [
                ;print "error:force out of bounds (high)"
                ;print xsum
                ;print i
                set xsum (2 * xsum / abs xsum)
          ]
        if abs ysum > 2
              [
                ;print "error:force out of bounds (high)"
                ;print ysum
                ;print i
                set ysum (2 * ysum / abs ysum)
          ]
        set xvect xvect + xsum
        set yvect yvect + ysum
      ]
    ]
     set i i + 1
  ]
  set i 1

end

to make-links
  set global-linkdist-list []
  set coeff-var-list []
  set global-link-count []
  while [j <= maximum ]
  [
    let x2 [xcor] of myotubes with [tubeno = j]
    ;only myotubes with +2 nuclei
    if length x2 > 2
    [
      set linkdistlist []
      ask myotubes with [tubeno = j]
      [
        let n1 0
        let n2 0
        let n3 0
       ;find index number of closest nuclei (Netlogo inicludes self in this list)
        if count my-links = 0
          [
            set largelink 0
            let aa [who] of min-n-of 3 myotubes with [tubeno = j] [distance myself]
            ;link nearest nuclei and record distance between
            set n1 item 0 aa
            set n2 item 1 aa
            set n3 item 2 aa
            ask myotube n1 [create-link-with myotube n2]
            ask link n1 n2
            [
              set largelink link-length
              set largest? true
              set force-check? false
            ]
            ask myotube n2 [create-link-with myotube n3]
            ask link n2 n3
            [
              set force-check? false
              ifelse link-length > largelink
              [
                set largelink link-length
                set largest? true
                ask link n1 n2
                [set largest? false]
              ]
              [set largest? false]
            ]
            ask myotube n1 [create-link-with myotube n3]
            ask link n1 n3
            [
              set force-check? false
              ifelse link-length > largelink
              [
                set largelink link-length
                set largest? true
                ask link n2 n3
                [set largest? false]
                ask link n1 n2
                [set largest? false]
              ]
              [set largest? false]
            ]
            ask link n1 n2
            [
              if largest? = false and force-check? = false
              [
                set force-check? true
                if GUI
                [set color green
                  set thickness radius / 2]
                set linkdistlist lput link-length linkdistlist
              ]

            ]
            ask link n2 n3
            [
              if largest? = false and force-check? = false
              [
                set force-check? true
                 if GUI
                [set color green
                  set thickness radius / 2]
                set linkdistlist lput link-length linkdistlist
              ]
            ]
            ask link n1 n3
            [
              if largest? = false and force-check? = false
              [
                set force-check? true
                if GUI
                [set color green
                  set thickness radius / 2]
                set linkdistlist lput link-length linkdistlist
              ]
            ]
            ask links with [largest? = true] [ die ]
        ]
      ]
      data-collect-local
    ]
    set j j + 1
  ]
  set j 1
end

to ptmove
  ask turtles
  [
    if tubeno = i
    [
    set xnew x_2 - xcent
    set ynew y_2 - ycent
    ;calculate points on linear regression line perpendicular to nucleus
    ;Calculate y-intercept (d) of perpendicular line
    ifelse m != 0
      [set d (ynew + (xnew / m))
    ;Calculate x-coordinate of point
        set xa (m * d) / (m * m + 1)
    ;Calculate y-coordinate of point
        set ya (m * xa)]
      [set xa xnew
        set ya 0]
      ;set lat-force as a random value from a given range
      let lat-force lat-force-min + random-float lat-force-diff

      let xRep-force ((xa - xnew) * lat-force)

      let yRep-force ((ya - ynew) * lat-force)
      if abs xRep-force > 1.5
      [
        ;print xRep-force
        set xRep-force (1.5 * xRep-force / abs xRep-force)
      ]
      if abs yRep-force > 1.5
      [
        ;print yRep-force
        set yRep-force (1.5 * yRep-force / abs yRep-force)
      ]
      set xvect xvect + xRep-force
      set yvect yvect + yRep-force

      set xRep-force 0
      set yRep-force 0
    ]
  ]
end
to showresticks
  ifelse show-resticks
  [
    ask turtles
    [
      set label tubeno
      ;set label fusion-age
      ;set label res-tick
      ;set label res-tot
      ;set label who
    ]
  ]
  [
    ask turtles
    [
    ;set label res-tick
      set label " "
    ]
  ]
end
to one-at-a-time
  if tubeset? = true
  [set tube-reset tubeno]
end

to myomove
  ;let y_wrap (ycor + yvect) mod max-pycor
  set y_2 (y_2 + yvect)
  set ycor (y_2 mod max-pycor)
  ;pendown
  ;let x_wrap (xcor + xvect) mod max-pxcor
  set x_2 (x_2 + xvect)
  set xcor (x_2 mod max-pxcor)

  if abs xvect > 1.5 and GUI
  [
    print who
    print "xvect"
    print xvect
  ]
   if abs yvect > 1.5 and GUI
  [
    print who
    print "yvect"
    print yvect
  ]
  set xvect 0
  set yvect 0

end

to data-collect-local
  let meanlinkdist mean linkdistlist
  ;show meanlinkdist
  set global-linkdist-list lput meanlinkdist global-linkdist-list
  let link-count (length linkdistlist + 1)
  set global-link-count lput link-count global-link-count
  if link-count > 2
  [
    let sdlinkdist standard-deviation linkdistlist
    ;show sdlinkdist
    ifelse meanlinkdist != 0
    [set coeff-var sdlinkdist / meanlinkdist]
    [set coeff-var 0]
    ;show coeff-var
    set coeff-var-list lput coeff-var coeff-var-list
    set sdlinkdist 0
    set coeff-var 0
  ]
  set meanlinkdist 0
  set link-count 0
end

to data-collect-global
  ;show length global-linkdist-list
  ;show mean global-link-count
  ;show standard-deviation global-link-count
  ;show mean global-linkdist-list
  ;show standard-deviation global-linkdist-list
  ;show median coeff-var-list
end
to-report number-myocytes
  report count myocytes
end
to-report number-myotubes
  report count myotubes ;length global-linkdist-list
end
to-report mean-nuc-per-tube
  report mean global-link-count
end
to-report sd-nuc-per-tube
  report standard-deviation global-link-count
end
to-report median-dist-nuc
  report median global-linkdist-list
end
to-report sd-dist-nuc
  report standard-deviation global-linkdist-list
end
to-report mean-coeff-var
  report mean coeff-var-list
end
@#$#@#$#@
GRAPHICS-WINDOW
0
0
1509
1510
-1
-1
1.5
1
20
1
1
1
0
1
1
1
0
1000
0
1000
1
1
1
ticks
30.0

BUTTON
32
14
99
47
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
22
114
118
147
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
1514
199
1686
232
myoblast-seed
myoblast-seed
0
5000
420.0
10
1
NIL
HORIZONTAL

SLIDER
1517
238
1689
271
angle
angle
0
90
0.0
1
1
NIL
HORIZONTAL

SLIDER
1704
239
1876
272
anglesd
anglesd
0
90
17.6
1
1
NIL
HORIZONTAL

SLIDER
1513
277
1685
310
speed
speed
0
10
1.25
0.01
1
NIL
HORIZONTAL

SLIDER
1700
279
1872
312
speedsd
speedsd
0
3
0.51
0.1
1
NIL
HORIZONTAL

SLIDER
1699
319
1872
352
radius
radius
0
60
20.0
1
1
NIL
HORIZONTAL

SLIDER
1515
361
1687
394
resmax
resmax
0
500
73.8
1
1
NIL
HORIZONTAL

PLOT
1507
524
1927
848
Cell count
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"myotube nuc" 1.0 0 -2674135 true "" "plotxy ticks count myotubes"
"myoblasts" 1.0 0 -13345367 true "" "plotxy ticks count myocytes"
"total cells" 1.0 0 -16777216 true "" "plot count turtles"
"fusion events" 1.0 0 -7500403 true "" "plotxy ticks fusecount"

SLIDER
1516
317
1688
350
Chance-division
Chance-division
0
1
7.51E-4
0.001
1
NIL
HORIZONTAL

PLOT
1508
854
1926
1225
Fusion Index
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plotxy ticks count myotubes / count turtles"

BUTTON
37
164
101
198
Rec
Rec
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
1703
200
1875
233
Myotube-seed
Myotube-seed
0
100
53.0
1
1
NIL
HORIZONTAL

SLIDER
1823
86
1944
119
nuc-forceB
nuc-forceB
0
20
0.0064
0.01
1
NIL
HORIZONTAL

SLIDER
1511
82
1646
115
lat-force-min
lat-force-min
0
100
0.02
1
1
NIL
HORIZONTAL

SWITCH
6
62
137
95
show-resticks
show-resticks
1
1
-1000

SWITCH
1804
37
1940
70
single-fusion
single-fusion
0
1
-1000

SWITCH
1634
38
1794
71
Nuclei_force
Nuclei_force
0
1
-1000

SLIDER
1509
138
1646
171
Max-MT-age
Max-MT-age
0
5
3.46
0.5
1
NIL
HORIZONTAL

SWITCH
1954
36
2057
69
GUI
GUI
1
1
-1000

SWITCH
1507
36
1620
69
DistDep
DistDep
1
1
-1000

SLIDER
1672
84
1808
117
lat-force-diff
lat-force-diff
0
100
0.07
1
1
NIL
HORIZONTAL

SLIDER
1701
361
1873
394
resmaxsd
resmaxsd
0
100
1.0
1
1
NIL
HORIZONTAL

SLIDER
1955
87
2077
120
nuc-forcesd
nuc-forcesd
0
0.01
0.0017
0.001
1
NIL
HORIZONTAL

SLIDER
1712
140
1884
173
Max-MT-age-SD
Max-MT-age-SD
0
5
1.15
1
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

This code example is a demo of a basic random walk.  At each step, the yellow turtle changes its heading randomly.

## THINGS TO NOTICE

The turtle's pen is down, so it leaves a trail behind it in the drawing.

## RELATED MODELS

Random Grid Walk Example: the same except that the random walk is constrained to lie on a grid

<!-- 2004 -->
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

myoblast
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 129 129 42

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
setup
repeat 2500 [ go ]
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10"/>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="anglesd">
      <value value="26"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="myoblast-seed">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resmax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="angle">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Chance-division">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="speed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="diameter">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="speedsd">
      <value value="0.35"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
