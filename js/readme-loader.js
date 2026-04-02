/*
 * readme-loader.js
 * Shared README renderer for CANS Lab project pages.
 * Hosted at: https://canslab1.github.io/js/readme-loader.js
 *
 * Auto-detects the repository name from the URL path.
 * Fetches ./README.md, renders it with marked.js, and highlights code blocks.
 */
(function () {
    'use strict';

    // Derive repo name from URL path: /RepoName/ → RepoName
    var pathParts = window.location.pathname.replace(/\/+$/, '').split('/');
    var repoName = pathParts[pathParts.length - 1] || 'unknown';

    // Sanitize renderer: strip dangerous tags (script, iframe, object, embed, form)
    var renderer = new marked.Renderer();
    var originalHtml = renderer.html;
    renderer.html = function (html) {
        var stripped = (typeof html === 'object' ? html.text || html.raw || '' : html)
            .replace(/<\s*\/?\s*(script|iframe|object|embed|form|link|meta|base)[^>]*>/gi, '');
        return stripped;
    };
    marked.use({ gfm: true, breaks: false, renderer: renderer });

    var container = document.getElementById('content');

    fetch('./README.md')
        .then(function (res) {
            if (!res.ok) throw new Error('HTTP ' + res.status);
            return res.text();
        })
        .then(function (md) {
            container.innerHTML = marked.parse(md);
            container.querySelectorAll('pre code').forEach(function (block) {
                hljs.highlightElement(block);
            });
        })
        .catch(function (err) {
            container.innerHTML =
                '<div class="error">' +
                '<p><strong>Could not load README.md</strong></p>' +
                '<p>' + err.message + '</p>' +
                '<p><a href="https://github.com/canslab1/' + repoName + '">View on GitHub</a></p>' +
                '</div>';
        });
})();
