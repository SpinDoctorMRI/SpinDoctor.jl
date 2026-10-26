import{_ as r,C as E,ap as g,o as l,c as e,j as i,a,ao as t,E as k,w as p,a2 as d}from"./chunks/framework.D9zB1E-I.js";const x=JSON.parse('{"title":"Apparent diffusion coefficients","description":"","frontmatter":{},"headers":[],"relativePath":"examples/generated/compare_adcs.md","filePath":"examples/generated/compare_adcs.md","lastUpdated":null}'),o={name:"examples/generated/compare_adcs.md"},y={class:"MathJax",jax:"SVG",overflow:"linebreak"},C={style:{"vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"4.965ex",height:"2.262ex",role:"img",focusable:"false",viewBox:"0 -750 2194.6 1000"},c={style:{"vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"9.505ex",height:"2.532ex",role:"img",focusable:"false",viewBox:"0 -869.3 4201.3 1119.3"},F={class:"vp-raw-html",innerHTML:`<div><div data-jscall-id="root" id="a65f9bb7-76fe-49e1-9d6a-a78ba6d5ae5c" class="bonito-fragment"><div><script src="bonito/js/Bonito.bundled10261065847571217372.js" type="module"><\/script><style></style></div><div><script type="module">Bonito.init_session("a65f9bb7-76fe-49e1-9d6a-a78ba6d5ae5c", Bonito.fetch_binary('bonito/bin/a278d6aa7b62a89ecba1ecc900c851216442cf5b-9978075562120054689.bin'), "root", false);
<\/script><div></div></div></div><div data-jscall-id="subsession-application-dom" id="f044797c-6142-4f22-9809-4115d87ccd60" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
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
</style></div><div><script type="module">Bonito.init_session("f044797c-6142-4f22-9809-4115d87ccd60", Bonito.fetch_binary('bonito/bin/27119247e313d7fdfc5e9bff9e24fce8220ce813-5286493482141051158.bin'), "sub", false);
<\/script><div data-jscall-id="1" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="2" height="500px" data-lm-suppress-shortcuts="true" width="600px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="3" class="wglmakie-spinner"></div></div></div></div></div>`},m={class:"vp-raw-html",innerHTML:`<div data-jscall-id="subsession-application-dom" id="dd3a5fbd-32c2-4e02-b566-f514117d8a3f" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
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
</style></div><div><script type="module">Bonito.init_session("dd3a5fbd-32c2-4e02-b566-f514117d8a3f", Bonito.fetch_binary('bonito/bin/dfe3e9e6b79d2e5a644d65ab555fd9d7201b6f28-6926893948567880894.bin'), "sub", false);
<\/script><div data-jscall-id="4" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="5" height="400px" data-lm-suppress-shortcuts="true" width="900px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="6" class="wglmakie-spinner"></div></div></div></div>`},u={class:"vp-raw-html",innerHTML:`<div data-jscall-id="subsession-application-dom" id="99338c4d-baf1-4b90-914b-67c472c99593" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
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
</style></div><div><script type="module">Bonito.init_session("99338c4d-baf1-4b90-914b-67c472c99593", Bonito.fetch_binary('bonito/bin/aeb9fce2e311d3f3d83303540bec005149b535fe-15398585789522286921.bin'), "sub", false);
<\/script><div data-jscall-id="7" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="8" height="350px" data-lm-suppress-shortcuts="true" width="650px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="9" class="wglmakie-spinner"></div></div></div></div>`};function b(f,s,D,v,A,B){const n=E("ClientOnly"),h=g("exec-scripts");return l(),e("div",null,[s[7]||(s[7]=i("h1",{id:"adc-methods",tabindex:"-1"},[a("Apparent diffusion coefficients "),i("a",{class:"header-anchor",href:"#adc-methods","aria-label":'Permalink to "Apparent diffusion coefficients {#adc-methods}"'},"​")],-1)),i("p",null,[s[3]||(s[3]=a("At low b-values, the signal attenuation is described by the apparent diffusion coefficient (ADC), ",-1)),i("mjx-container",y,[(l(),e("svg",C,[...s[0]||(s[0]=[t("",1)])])),s[2]||(s[2]=i("mjx-break",{size:"4"}," ",-1)),(l(),e("svg",c,[...s[1]||(s[1]=[t("",1)])]))]),s[4]||(s[4]=a(" (see ",-1)),s[5]||(s[5]=i("a",{href:"/SpinDoctor.jl/previews/PR48/theory/adc#adc"},"Apparent diffusion coefficient",-1)),s[6]||(s[6]=a("). SpinDoctor computes the ADC in several ways:",-1))]),s[8]||(s[8]=t("",6)),k(n,null,{default:p(()=>[d(i("div",F,null,512),[[h]])]),_:1}),s[9]||(s[9]=t("",21)),k(n,null,{default:p(()=>[d(i("div",m,null,512),[[h]])]),_:1}),s[10]||(s[10]=t("",4)),k(n,null,{default:p(()=>[d(i("div",u,null,512),[[h]])]),_:1}),s[11]||(s[11]=i("hr",null,null,-1)),s[12]||(s[12]=i("p",null,[i("em",null,[a("This page was generated using "),i("a",{href:"https://github.com/fredrikekre/Literate.jl",target:"_blank",rel:"noreferrer"},"Literate.jl"),a(".")])],-1))])}const w=r(o,[["render",b]]);export{x as __pageData,w as default};
