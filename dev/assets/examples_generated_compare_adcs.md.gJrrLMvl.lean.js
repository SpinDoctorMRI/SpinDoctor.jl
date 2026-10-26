import{_ as r,C as E,ap as g,o as l,c as e,j as i,a,ao as t,E as k,w as p,a2 as d}from"./chunks/framework.CXYRQASs.js";const x=JSON.parse('{"title":"Apparent diffusion coefficients","description":"","frontmatter":{},"headers":[],"relativePath":"examples/generated/compare_adcs.md","filePath":"examples/generated/compare_adcs.md","lastUpdated":null}'),o={name:"examples/generated/compare_adcs.md"},y={class:"MathJax",jax:"SVG",overflow:"linebreak"},C={style:{"vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"4.965ex",height:"2.262ex",role:"img",focusable:"false",viewBox:"0 -750 2194.6 1000"},c={style:{"vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"9.505ex",height:"2.532ex",role:"img",focusable:"false",viewBox:"0 -869.3 4201.3 1119.3"},F={class:"vp-raw-html",innerHTML:`<div><div data-jscall-id="root" id="dd041f95-72dd-4f88-9fce-8f686003266b" class="bonito-fragment"><div><script src="bonito/js/Bonito.bundled10261065847571217372.js" type="module"><\/script><style></style></div><div><script type="module">Bonito.init_session("dd041f95-72dd-4f88-9fce-8f686003266b", Bonito.fetch_binary('bonito/bin/1def173572db8d90b692bbd17f02c547dd9a812d-9051008517232164107.bin'), "root", false);
<\/script><div></div></div></div><div data-jscall-id="subsession-application-dom" id="47ead084-f094-4012-97dc-b98cab251d14" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
  animation: wglmakie-spin 1.0s linear infinite;
  position: absolute;
  left: 50%;
  z-index: 1000;
  border-radius: 50%;
  width: 40px;
  border: 4px solid rgba(0, 0, 0, 0.1);
  transform: translate(-50%, -50%);
  border-top: 4px solid currentColor;
  height: 40px;
  pointer-events: none;
  top: 50%;
}
@keyframes wglmakie-spin {
  100% {
    transform: translate(-50%, -50%) rotate(360deg);
  }
  0% {
    transform: translate(-50%, -50%) rotate(0deg);
  }
}
</style></div><div><script type="module">Bonito.init_session("47ead084-f094-4012-97dc-b98cab251d14", Bonito.fetch_binary('bonito/bin/8d4fe09dff6df0f57c73e3813718302dc84f88f2-4868992201035942710.bin'), "sub", false);
<\/script><div data-jscall-id="1" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="2" height="500px" data-lm-suppress-shortcuts="true" width="600px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="3" class="wglmakie-spinner"></div></div></div></div></div>`},m={class:"vp-raw-html",innerHTML:`<div data-jscall-id="subsession-application-dom" id="453050e8-5010-48b7-be21-2940e4a51c84" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
  animation: wglmakie-spin 1.0s linear infinite;
  position: absolute;
  left: 50%;
  z-index: 1000;
  border-radius: 50%;
  width: 40px;
  border: 4px solid rgba(0, 0, 0, 0.1);
  transform: translate(-50%, -50%);
  border-top: 4px solid currentColor;
  height: 40px;
  pointer-events: none;
  top: 50%;
}
@keyframes wglmakie-spin {
  100% {
    transform: translate(-50%, -50%) rotate(360deg);
  }
  0% {
    transform: translate(-50%, -50%) rotate(0deg);
  }
}
</style></div><div><script type="module">Bonito.init_session("453050e8-5010-48b7-be21-2940e4a51c84", Bonito.fetch_binary('bonito/bin/58ea4081f643bcad048f734af370c6a4c5ed6aa6-16821067196446004697.bin'), "sub", false);
<\/script><div data-jscall-id="4" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="5" height="400px" data-lm-suppress-shortcuts="true" width="900px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="6" class="wglmakie-spinner"></div></div></div></div>`},u={class:"vp-raw-html",innerHTML:`<div data-jscall-id="subsession-application-dom" id="89148d3e-d516-4b0d-99d7-c72d6d5450e0" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
  animation: wglmakie-spin 1.0s linear infinite;
  position: absolute;
  left: 50%;
  z-index: 1000;
  border-radius: 50%;
  width: 40px;
  border: 4px solid rgba(0, 0, 0, 0.1);
  transform: translate(-50%, -50%);
  border-top: 4px solid currentColor;
  height: 40px;
  pointer-events: none;
  top: 50%;
}
@keyframes wglmakie-spin {
  100% {
    transform: translate(-50%, -50%) rotate(360deg);
  }
  0% {
    transform: translate(-50%, -50%) rotate(0deg);
  }
}
</style></div><div><script type="module">Bonito.init_session("89148d3e-d516-4b0d-99d7-c72d6d5450e0", Bonito.fetch_binary('bonito/bin/13fd99e3388ff59f08dc5fe7822943ef8ba46ad7-17230820295119550802.bin'), "sub", false);
<\/script><div data-jscall-id="7" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="8" height="350px" data-lm-suppress-shortcuts="true" width="650px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="9" class="wglmakie-spinner"></div></div></div></div>`};function b(f,s,D,v,A,B){const n=E("ClientOnly"),h=g("exec-scripts");return l(),e("div",null,[s[7]||(s[7]=i("h1",{id:"adc-methods",tabindex:"-1"},[a("Apparent diffusion coefficients "),i("a",{class:"header-anchor",href:"#adc-methods","aria-label":'Permalink to "Apparent diffusion coefficients {#adc-methods}"'},"​")],-1)),i("p",null,[s[3]||(s[3]=a("At low b-values, the signal attenuation is described by the apparent diffusion coefficient (ADC), ",-1)),i("mjx-container",y,[(l(),e("svg",C,[...s[0]||(s[0]=[t("",1)])])),s[2]||(s[2]=i("mjx-break",{size:"4"}," ",-1)),(l(),e("svg",c,[...s[1]||(s[1]=[t("",1)])]))]),s[4]||(s[4]=a(" (see ",-1)),s[5]||(s[5]=i("a",{href:"/SpinDoctor.jl/dev/theory/adc#adc"},"Apparent diffusion coefficient",-1)),s[6]||(s[6]=a("). SpinDoctor computes the ADC in several ways:",-1))]),s[8]||(s[8]=t("",6)),k(n,null,{default:p(()=>[d(i("div",F,null,512),[[h]])]),_:1}),s[9]||(s[9]=t("",21)),k(n,null,{default:p(()=>[d(i("div",m,null,512),[[h]])]),_:1}),s[10]||(s[10]=t("",4)),k(n,null,{default:p(()=>[d(i("div",u,null,512),[[h]])]),_:1}),s[11]||(s[11]=i("hr",null,null,-1)),s[12]||(s[12]=i("p",null,[i("em",null,[a("This page was generated using "),i("a",{href:"https://github.com/fredrikekre/Literate.jl",target:"_blank",rel:"noreferrer"},"Literate.jl"),a(".")])],-1))])}const w=r(o,[["render",b]]);export{x as __pageData,w as default};
