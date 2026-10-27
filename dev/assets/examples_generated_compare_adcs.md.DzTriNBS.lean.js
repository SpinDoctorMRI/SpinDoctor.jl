import{_ as r,C as E,ap as g,o as l,c as e,j as i,a,ao as t,E as k,w as p,a2 as d}from"./chunks/framework.CXYRQASs.js";const x=JSON.parse('{"title":"Apparent diffusion coefficients","description":"","frontmatter":{},"headers":[],"relativePath":"examples/generated/compare_adcs.md","filePath":"examples/generated/compare_adcs.md","lastUpdated":null}'),o={name:"examples/generated/compare_adcs.md"},y={class:"MathJax",jax:"SVG",overflow:"linebreak"},C={style:{"vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"4.965ex",height:"2.262ex",role:"img",focusable:"false",viewBox:"0 -750 2194.6 1000"},c={style:{"vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"9.505ex",height:"2.532ex",role:"img",focusable:"false",viewBox:"0 -869.3 4201.3 1119.3"},F={class:"vp-raw-html",innerHTML:`<div><div data-jscall-id="root" id="da7cd5f4-abba-4ff4-b46c-0b737f23e118" class="bonito-fragment"><div><script src="bonito/js/Bonito.bundled10261065847571217372.js" type="module"><\/script><style></style></div><div><script type="module">Bonito.init_session("da7cd5f4-abba-4ff4-b46c-0b737f23e118", Bonito.fetch_binary('bonito/bin/b41b5e2dac90b646e5e8ebb429a57584253cc485-10752713237312383751.bin'), "root", false);
<\/script><div></div></div></div><div data-jscall-id="subsession-application-dom" id="787ea968-155e-4047-818a-4c2a4d0d5b99" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
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
</style></div><div><script type="module">Bonito.init_session("787ea968-155e-4047-818a-4c2a4d0d5b99", Bonito.fetch_binary('bonito/bin/d5dd57f0f7b00ffa3005b97828f1c1eed09610f7-5745042079962044837.bin'), "sub", false);
<\/script><div data-jscall-id="1" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="2" height="500px" data-lm-suppress-shortcuts="true" width="600px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="3" class="wglmakie-spinner"></div></div></div></div></div>`},m={class:"vp-raw-html",innerHTML:`<div data-jscall-id="subsession-application-dom" id="5b2fedea-2e71-4ff6-9839-bc9d24efd5f5" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
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
</style></div><div><script type="module">Bonito.init_session("5b2fedea-2e71-4ff6-9839-bc9d24efd5f5", Bonito.fetch_binary('bonito/bin/5e969dbf00c554a6084f6253d4f2431a37d87d07-11107403669057771381.bin'), "sub", false);
<\/script><div data-jscall-id="4" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="5" height="400px" data-lm-suppress-shortcuts="true" width="900px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="6" class="wglmakie-spinner"></div></div></div></div>`},u={class:"vp-raw-html",innerHTML:`<div data-jscall-id="subsession-application-dom" id="3c5ab8b4-043f-4845-9d43-ff9766462b9f" class="bonito-fragment"><div><style></style><style>.wglmakie-spinner {
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
</style></div><div><script type="module">Bonito.init_session("3c5ab8b4-043f-4845-9d43-ff9766462b9f", Bonito.fetch_binary('bonito/bin/df6f4cd4e1df2742e2cad0b72f5cddab2daa25d5-8074352285205372841.bin'), "sub", false);
<\/script><div data-jscall-id="7" style="width: 100%; height: 100%; position: relative;"><canvas data-jscall-id="8" height="350px" data-lm-suppress-shortcuts="true" width="650px" data-jp-suppress-context-menu tabindex="0" style="display: block"></canvas><div data-jscall-id="9" class="wglmakie-spinner"></div></div></div></div>`};function b(f,s,D,v,A,B){const n=E("ClientOnly"),h=g("exec-scripts");return l(),e("div",null,[s[7]||(s[7]=i("h1",{id:"adc-methods",tabindex:"-1"},[a("Apparent diffusion coefficients "),i("a",{class:"header-anchor",href:"#adc-methods","aria-label":'Permalink to "Apparent diffusion coefficients {#adc-methods}"'},"​")],-1)),i("p",null,[s[3]||(s[3]=a("At low b-values, the signal attenuation is described by the apparent diffusion coefficient (ADC), ",-1)),i("mjx-container",y,[(l(),e("svg",C,[...s[0]||(s[0]=[t("",1)])])),s[2]||(s[2]=i("mjx-break",{size:"4"}," ",-1)),(l(),e("svg",c,[...s[1]||(s[1]=[t("",1)])]))]),s[4]||(s[4]=a(" (see ",-1)),s[5]||(s[5]=i("a",{href:"/SpinDoctor.jl/dev/theory/adc#adc"},"Apparent diffusion coefficient",-1)),s[6]||(s[6]=a("). SpinDoctor computes the ADC in several ways:",-1))]),s[8]||(s[8]=t("",6)),k(n,null,{default:p(()=>[d(i("div",F,null,512),[[h]])]),_:1}),s[9]||(s[9]=t("",21)),k(n,null,{default:p(()=>[d(i("div",m,null,512),[[h]])]),_:1}),s[10]||(s[10]=t("",4)),k(n,null,{default:p(()=>[d(i("div",u,null,512),[[h]])]),_:1}),s[11]||(s[11]=i("hr",null,null,-1)),s[12]||(s[12]=i("p",null,[i("em",null,[a("This page was generated using "),i("a",{href:"https://github.com/fredrikekre/Literate.jl",target:"_blank",rel:"noreferrer"},"Literate.jl"),a(".")])],-1))])}const w=r(o,[["render",b]]);export{x as __pageData,w as default};
