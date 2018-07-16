/* Interactive Simulated Sonification of the Double Pendulum */

/* Run the following two synth defs BEFORE executing the main program loop. */
(
SynthDef.new(\custom, {
	arg freq=263, pan=0, mod=0, amp=1, mul=1;
	var sig;
	sig = Pan2.ar(SinOsc.ar(freq+SinOsc.kr(mod, mul:amp), mul:0.25)*mul, pan);
	Out.ar(0, sig);
}).add;
)

(
SynthDef.new(\custperc, {
	arg freq=263, on=0, pan=0, nHarm = 1, mul=1;
	var sig;
	sig = Pan2.ar(Blip.ar(freq, nHarm, mul:on)*mul, pan);
	Out.ar(0, sig);
}).add;
)


/* The following code is the main program loop. */
(
var w, b1, gslider, grange, gtext, gvalue, tslider, trange, ttext, tvalue, l1slider, l1range, l1text, l1value, l2slider, l2range, l2text, l2value, m1slider, m1range, m1text, m1value, m2slider, m2range, m2text, m2value, run = true;

/* Simulation Parameters */

~g = 9.81; //gravity in m/s^2
~l1 = 1; //length of pendulum 1 in m
~l2 = 1;//length of pendulum 2 in m
~m1 = ~l1; //mass of pendulum 1 in kg
~m2 = ~l2; //mass of pendulum 2 in kg
~th1 = pi.rand; //initial angle in radians
~w1 = 0;  //initial angular speed in radians/s
~th2 = pi.rand; //initial angle in radians
~w2 = 0;  //initial angular speed in radians/s
~h = 0.001; //stepsize of integration
~t = 9; //speed of simulation (# steps per update)

/*
Simulation Routine (everytime it is called, the internal loop runs exactly one time)
*/

p = Routine({
	var x, y, yin, del, den1, den2, f, k1, k2, k3, k4;

	yin = [~th1, ~w1, ~th2, ~w2];
	f = {
		arg yin, del, den1, den2;
		[yin[1], (((~m2*~l1*yin[1]*yin[1]*(del.sin)*(del.cos))		+ (~m2*~g*sin(yin[2])*(del.cos)) + (~m2*~l1*yin[3]*yin[3]*(del.sin))	- ((~m1+~m1)*~g*(yin[0].sin)))/den1), (yin[3]), (((-1*~m2*~l2*yin[3]*yin[3]*(del.sin)*cos(del.cos))+ ((~m1+~m2)*~g*(yin[0].sin)*(del.cos))- ((~m1+~m2)*~l1*yin[1]*yin[1]*(del.sin))- ((~m1+~m2)*~g*(yin[2].sin)))/den2)]
	};

	//fourth order Runge-Kutta (dy/dx = f)
	loop({
		yin = [~th1, ~w1, ~th2, ~w2];

		/*
		The nested "for" loop allows simulation to advance
		more (i.e. run faster) for larger ~t (advances by
		~h*~t steps every ~h*10 seconds). Values larger
		than about 100-500 freeze up the server and require
		a reboot of the interpreter, so be careful/
		*/

		for(1, ~t,{
			del = yin[2]-yin[0];
			den1 = (~m1+~m2)*~l1 - (~m2*~l1*(del.cos)*(del.cos));
			den2 = (~l2/~l1)*den1;

			//first step

			k1 = ~h*(f.valueArray(yin, del, den1, den2));

			//second step

			yin = yin + (0.5*k1);
			del = yin[2]-yin[0];
			den1 = (~m1+~m2)*~l1 - (~m2*~l1*(del.cos)*(del.cos));
			den2 = (~l2/~l1)*den1;
			k2 = ~h*(f.valueArray(yin, del, den1, den2));

			//third step

			yin = yin + (0.5*k2) - (0.5*k1);
			del = yin[2]-yin[0];
			den1 = (~m1+~m2)*~l1 - (~m2*~l1*(del.cos)*(del.cos));
			den2 = (~l2/~l1)*den1;
			k3 = ~h*(f.valueArray(yin, del, den1, den2));

			//fourth step

			yin = yin + k3 - (0.5*k2);
			del = yin[2]-yin[0];
			den1 = (~m1+~m2)*~l1 - (~m2*~l1*(del.cos)*(del.cos));
			den2 = (~l2/~l1)*den1;
			k4 = ~h*(f.valueArray(yin, del, den1, den2));
			yin = yin - k3;

			//final value of y for the nth step

			yin = yin + (k1/6) + (k2/3) + (k3/3) + (k4/6);
			yin[0] = yin[0]%(2*pi);//range between 0 and 2pi
			yin[2] = yin[2]%(2*pi);//range between 0 and 2pi

		});

		/* Pendulum has now advanced ~h*~t steps */

		~th1 = yin[0];
		~w1 = yin[1];
		~th2 = yin[2];
		~w2 = yin[3];
		//yin = [~th1, ~w1, ~th2, ~w2];

		//conversion of polar to cartesian coordinates

		x = ~l1*(~th1.sin) + (~l2*((~th2).sin));
		y = (~l1+~l2)-(1*~l1*(~th1.cos)) - (~l2*((~th2).cos));

		//map of phase space coordinates to sound parameters

		~freq1 = 112.5 *( 2**(y));
		~pan1 = x;
		~mod1 = ((~w2.abs)+(~w1.abs))*5*(~th1.sin);
		~amp1 = ~freq1/10;
		~freq2 = 112.5 *( 3**(y));
		~pan2 = x;
		~mod2 = ((~w2.abs)+(~w1.abs))*5*(~th1.sin);
		~amp2 = ~freq2/10;
		~freq3 = ~freq1;
		~pan3 = x;
		~nH3 = 1+(20*((~th1+~th2).sin));
		if (~th1<0.1, {~on3 = 1}, {~on3 = 0});
		~freq4 = ~freq2;
		~pan4 = x;
		~nH4 = 1+(20*(~th2.sin));
		if ((~th1-~th2).abs<0.1, {~on4 = 1}, {~on4 = 0});
		[y, x].yield;
	});
});

/* Drawing Function */

w = Window.new(bounds:Rect(400,300,550,550)).front;
w.drawFunc = {
	Pen.use{
		var r = Rect.newSides(-5,0,5,100*~l1);
		Pen.fillColor = Color.blue;
		Pen.translate(275,210);
		Pen.rotate(-1*~th1);
		Pen.fillRect(r);
		r = Rect.newSides(-5,0,5,100*~l2);
		Pen.fillColor = Color.red;
		Pen.translate(0,100*~l1);
		Pen.rotate(~th1);
		Pen.rotate(-1*~th2);
		Pen.fillRect(r);
		Pen.rotate(~th2);
		}
};

/* Control Panel */

//gravity slider

grange = ControlSpec(0.1, 20, \linear, 0.01);
gtext = StaticText(w, Rect(75,415,50,30));
gtext.string = "Grav. (m/s^2)";
gvalue = NumberBox(w, Rect(20,420,45,24));
gslider = Slider(w, Rect(20,450,114,24)).value_(0.4873).action_({
	gvalue.value_(grange.map(gslider.value));
	~g = gvalue.value;
});
gslider.action.value;
//speed slider
trange = ControlSpec(1, 90, \exponential, 0.5);
ttext = StaticText(w, Rect(75,480,50,30));
ttext.string = "Sim. Speed";
tvalue = NumberBox(w, Rect(20,480,45,24));
tslider = Slider(w, Rect(20,510,114,24)).value_(0.1).action_({
	tvalue.value_(trange.map(tslider.value));
	~t = tvalue.value;
});
//l1 slider
l1range = ControlSpec(0.05, 1, \linear, 0.01);
l1text = StaticText(w, Rect(220,415,50,30));
l1text.string = "L1 (m)";
l1value = NumberBox(w, Rect(165,420,45,24));
l1slider = Slider(w, Rect(165,450,114,24)).value_(1).action_({
	l1value.value_(l1range.map(l1slider.value));
	~l1 = l1value.value;
});
//l2slider
l2range = ControlSpec(0.05, 1, \linear, 0.01);
l2text = StaticText(w, Rect(220,480,50,30));
l2text.string = "L2 (m)";
l2value = NumberBox(w, Rect(165,480,45,24));
l2slider = Slider(w, Rect(165,510,114,24)).value_(1).action_({
	l2value.value_(l2range.map(l2slider.value));
	~l2 = l2value.value;
});
//m1slider
m1range = ControlSpec(0.05, 1, \linear, 0.01);
m1text = StaticText(w, Rect(375,420,50,30));
m1text.string = "M1 (kg)";
m1value = NumberBox(w, Rect(310,420,45,24));
m1slider = Slider(w, Rect(310,450,114,24)).value_(1).action_({
	m1value.value_(m1range.map(m1slider.value));
	~m1 = m1value.value;
});
//m2slider
m2range = ControlSpec(0.05, 1, \linear, 0.01);
m2text = StaticText(w, Rect(375,480,50,30));
m2text.string = "M2 (kg)";
m2value = NumberBox(w, Rect(310,480,45,24));
m2slider = Slider(w, Rect(310,510,114,24)).value_(1).action_({
	m2value.value_(m2range.map(m2slider.value));
	~m2 = m2value.value;
});
//reset button
b1 = Button.new(w, Rect(450, 470, 69, 24));
b1.states_([["Reset"]]);
b1.action_({arg button;
	if(button.value == 0, {
		~th1 = (2*pi).rand;
		~w1 = 0;
		~th2 = (2*pi).rand;
		~w2 = 0;
	});
});
/* Sound Generation */

~z1 = Synth.new(\custom,[\mul, 0.5]);
~z2 = Synth.new(\custom,[\mul, 0.5]);
~z3 = Synth.new(\custperc,[\mul, 0.5]);
~z4 = Synth.new(\custperc,[\mul, 0.5]);
// s.scope;

/*
This while loop initiates a thread which continually updates the pendulum simulation (not graphics)
*/

{while{run}{
	~z1.set(
		\freq, ~freq1,
		\pan, ~pan1,
		\mod, ~mod1,
		\amp, ~amp1,
	);
	~z2.set(
		\freq, ~freq2,
		\pan, ~pan2,
		\mod, ~mod2,
		\amp, ~amp2,
	);
	~z3.set(
		\freq, ~freq3,
		\pan, ~pan3,
		\on, ~on3,
		\nHarm, ~nH3,
	);
	~z4.set(
		\freq, ~freq4,
		\pan, ~pan4,
		\on, ~on4,
		\nHarm, ~nH4,
	);
	p.value;
	(~h*10).wait;
}
}.fork(SystemClock);

/*
The choice of wait time (proportional to ~h) ensures that changes in ~h don't affect the speed of the simulation.
*/

/*
This while loop initiates a thread which continually updates the drawing of the pendulum
*/

{while{run}{
	w.refresh;
	(1/60).wait;
}
}.fork(AppClock);
)

(
x = Ndef(\osc, {
	|freq, pan|
	Out.ar(0, SinOsc.ar(freq));
}).add;

NdefGui(x);
)
