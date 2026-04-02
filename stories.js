/* ===== Microsoft Clarity ===== */
(function(c,l,a,r,i,t,y){
    c[a]=c[a]||function(){(c[a].q=c[a].q||[]).push(arguments)};
    t=l.createElement(r);t.async=1;t.src="https://www.clarity.ms/tag/"+i+"?ref=bwt";
    y=l.getElementsByTagName(r)[0];y.parentNode.insertBefore(t,y);
})(window, document, "clarity", "script", "rzlnthqbys");

/* ===== Mobile Menu Toggle ===== */
function toggleMenu() {
    const nav = document.getElementById("navLinks");
    if (nav) nav.classList.toggle("show");
}

/* ===== Table of Contents Toggle ===== */
function toggleTOC() {
    const tocContent = document.getElementById('tocContent');
    const arrow = document.querySelector('.toc-arrow');
    const hint = document.querySelector('.toc-hint');
    if (!tocContent || !arrow || !hint) return;

    const tocBtn = document.querySelector('.toc-title');
    const isEn = document.documentElement.lang === 'en';
    if (tocContent.classList.contains('expanded')) {
        tocContent.classList.remove('expanded');
        arrow.classList.remove('expanded');
        hint.textContent = isEn ? 'Tap to expand' : '點擊展開';
        if (tocBtn) tocBtn.setAttribute('aria-expanded', 'false');
    } else {
        tocContent.classList.add('expanded');
        arrow.classList.add('expanded');
        hint.textContent = isEn ? 'Tap to collapse' : '點擊收起';
        if (tocBtn) tocBtn.setAttribute('aria-expanded', 'true');
    }
}

/* ===== Back to Top ===== */
window.addEventListener('scroll', function () {
    const backToTop = document.querySelector('.back-to-top');
    if (!backToTop) return;
    if (window.pageYOffset > 300) {
        backToTop.classList.add('visible');
    } else {
        backToTop.classList.remove('visible');
    }
}, { passive: true });

function scrollToTop() {
    window.scrollTo({
        top: 0,
        behavior: 'smooth'
    });
}

/* ===== Smooth Anchor Scrolling ===== */
function smoothScrollToAnchor(anchorElement) {
    const targetId = anchorElement.getAttribute('href');
    const target = document.querySelector(targetId);
    if (target) {
        setTimeout(() => {
            target.scrollIntoView({ behavior: 'smooth', block: 'start' });
            history.pushState(null, '', targetId);
        }, 300);
    }
}

document.addEventListener('DOMContentLoaded', function () {
    /* Bind nav toggle button */
    const navToggle = document.getElementById('navToggle');
    if (navToggle) navToggle.addEventListener('click', toggleMenu);

    /* Bind TOC toggle button */
    const tocToggle = document.getElementById('tocToggle');
    if (tocToggle) tocToggle.addEventListener('click', toggleTOC);

    /* Bind back-to-top button */
    const backToTop = document.getElementById('backToTop');
    if (backToTop) backToTop.addEventListener('click', scrollToTop);

    /* Smooth anchor scrolling */
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function (e) {
            if (this.getAttribute('href').length > 1) {
                e.preventDefault();
                smoothScrollToAnchor(this);
            }
        });
    });
});

/* ===== Fade-in Animation Observer ===== */
document.addEventListener('DOMContentLoaded', function () {
    const observer = new IntersectionObserver((entries) => {
        entries.forEach((entry) => {
            if (entry.isIntersecting) {
                entry.target.style.animationDelay = '0s';
                entry.target.style.animationPlayState = 'running';
            }
        });
    });

    document.querySelectorAll('.chapter, .section').forEach((el) => {
        observer.observe(el);
    });
});
