function showTool(toolId, btn) {

    document.querySelectorAll('.tool').forEach(t => t.classList.remove('active'));

    document.getElementById(toolId).classList.add('active');

    document.querySelectorAll('.topnav button, .dropdown-content button').forEach(b => b.classList.remove('active'));

    btn.classList.add('active');
}

var PHYLODENDRON_CONFIG = window.PHYLODENDRON_CONFIG || {};
var newickViewer;
var newickRawTreeString = PHYLODENDRON_CONFIG.treeNewick || null;

window.enableSplitNewickLabelColors = function() {
    if (!window.Phylocanvas || !Phylocanvas.Branch || !Phylocanvas.Branch.prototype) {
        return;
    }

    var branchProto = Phylocanvas.Branch.prototype;
    if (branchProto.__splitLabelColorsPatched) {
        return;
    }

    var originalDrawBranchLabels = branchProto.drawBranchLabels;

    branchProto.drawBranchLabels = function() {
        if (!this.tree || !this.tree.customLabelColours) {
            return originalDrawBranchLabels.call(this);
        }

        this.canvas.save();
        var labelStyle = this.internalLabelStyle || this.tree.internalLabelStyle;
        this.canvas.font = labelStyle.format + ' ' + labelStyle.textSize + 'pt ' + labelStyle.font;
        this.canvas.textBaseline = 'middle';
        this.canvas.textAlign = 'center';
        var em = this.canvas.measureText('M').width * 2 / 3;

        var x = this.tree.type.branchScalingAxis === 'y' ? this.centerx : (this.startx + this.centerx) / 2;
        var y = this.tree.type.branchScalingAxis === 'x' ? this.centery : (this.starty + this.centery) / 2;

        if (this.tree.showBranchLengthLabels && this.tree.branchLengthLabelPredicate(this)) {
            this.canvas.fillStyle = this.tree.newickBranchLengthColour || labelStyle.colour;
            this.canvas.fillText(this.branchLength.toFixed(2), x, y + em);
        }

        if (this.tree.showInternalNodeLabels && !this.leaf && this.label) {
            this.canvas.fillStyle = this.tree.newickInternalNodeColour || labelStyle.colour;
            this.canvas.fillText(this.label, x, y - em);
        }

        this.canvas.restore();
    };

    branchProto.__splitLabelColorsPatched = true;
};

function getNewickStyleValues() {
    return {
        branchColor: (document.getElementById('newick-branch-color') || {}).value || '#111827',
        leafColor: (document.getElementById('newick-leaf-color') || {}).value || '#111827',
        branchLengthColor: (document.getElementById('newick-branch-length-color') || {}).value || '#111827',
        internalNodeColor: (document.getElementById('newick-internal-node-color') || {}).value || '#111827'
    };
}

function setTreeAnnotations(prefix) {
    var titleInput = document.getElementById(prefix + '-title-input');
    var xInput = document.getElementById(prefix + '-x-input');
    var yInput = document.getElementById(prefix + '-y-input');

    var titleLabel = document.getElementById(prefix + '-title-label');
    var xLabel = document.getElementById(prefix + '-x-label');
    var yLabel = document.getElementById(prefix + '-y-label');

    if (titleLabel) {
        titleLabel.textContent = titleInput ? titleInput.value.trim() : '';
    }
    if (xLabel) {
        xLabel.textContent = xInput ? xInput.value.trim() : '';
    }
    if (yLabel) {
        yLabel.textContent = yInput ? yInput.value.trim() : '';
    }
}

window.toggleControlCard = function(button) {
    var card = button ? button.closest('.newick-control-card') : null;
    if (!card) {
        return;
    }

    var shouldCollapse = !card.classList.contains('is-collapsed');
    card.classList.toggle('is-collapsed', shouldCollapse);

    var expanded = !shouldCollapse;
    button.setAttribute('aria-expanded', expanded ? 'true' : 'false');

    var icon = button.querySelector('.newick-card-toggle-icon');
    if (icon) {
        icon.textContent = expanded ? '−' : '+';
    }

    var workspace = card.closest('.newick-workspace');
    if (!workspace) {
        return;
    }

    // Wait one frame so the collapsed layout is applied before measuring height.
    window.requestAnimationFrame(function() {
        if (workspace.id === 'bayes-workspace') {
            syncInferenceWorkspaceHeights('bayes');
        } else if (workspace.id === 'iqtree-workspace') {
            syncInferenceWorkspaceHeights('iqtree');
        } else if (workspace.id === 'mpboot-workspace') {
            syncInferenceWorkspaceHeights('mpboot');
        } else if (workspace.id === 'distance-workspace') {
            syncInferenceWorkspaceHeights('distance');
        } else {
            syncNewickWorkspaceHeights();
        }
    });
};

window.updateNewickAnnotations = function() {
    setTreeAnnotations('newick');
};

function getTreeAnnotationValues(prefix) {
    var titleInput = document.getElementById(prefix + '-title-input');
    var xInput = document.getElementById(prefix + '-x-input');
    var yInput = document.getElementById(prefix + '-y-input');

    return {
        title: titleInput ? titleInput.value.trim() : '',
        xLabel: xInput ? xInput.value.trim() : '',
        yLabel: yInput ? yInput.value.trim() : ''
    };
}

function buildAnnotatedExportCanvas(sourceCanvas, annotationPrefix) {
    var annotations = getTreeAnnotationValues(annotationPrefix);
    var hasTitle = annotations.title.length > 0;
    var hasXLabel = annotations.xLabel.length > 0;
    var hasYLabel = annotations.yLabel.length > 0;

    var topPad = hasTitle ? 90 : 24;
    var bottomPad = hasXLabel ? 90 : 24;
    var leftPad = hasYLabel ? 120 : 24;
    var rightPad = 24;

    var outCanvas = document.createElement('canvas');
    outCanvas.width = sourceCanvas.width + leftPad + rightPad;
    outCanvas.height = sourceCanvas.height + topPad + bottomPad;

    var ctx = outCanvas.getContext('2d');
    ctx.fillStyle = '#FFFFFF';
    ctx.fillRect(0, 0, outCanvas.width, outCanvas.height);
    ctx.drawImage(sourceCanvas, leftPad, topPad);

    ctx.fillStyle = '#111827';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    if (hasTitle) {
        ctx.font = '700 48px "Segoe UI", Roboto, Arial, sans-serif';
        ctx.fillText(annotations.title, outCanvas.width / 2, topPad / 2);
    }

    if (hasXLabel) {
        ctx.font = '600 34px "Segoe UI", Roboto, Arial, sans-serif';
        ctx.fillText(annotations.xLabel, outCanvas.width / 2, outCanvas.height - (bottomPad / 2));
    }

    if (hasYLabel) {
        ctx.save();
        ctx.translate(leftPad / 2, outCanvas.height / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.font = '600 34px "Segoe UI", Roboto, Arial, sans-serif';
        ctx.fillText(annotations.yLabel, 0, 0);
        ctx.restore();
    }

    return outCanvas;
}

window.updateNewickStyles = function() {
    if (!newickViewer) {
        return;
    }

    var colors = getNewickStyleValues();

    newickViewer.branchColour = colors.branchColor;
    newickViewer.customLabelColours = true;
    newickViewer.newickBranchLengthColour = colors.branchLengthColor;
    newickViewer.newickInternalNodeColour = colors.internalNodeColor;
    newickViewer.internalLabelStyle = newickViewer.internalLabelStyle || {};
    newickViewer.internalLabelStyle.colour = colors.branchLengthColor;

    if (newickViewer.branches) {
        Object.keys(newickViewer.branches).forEach(function(branchId) {
            var branch = newickViewer.branches[branchId];
            if (!branch) {
                return;
            }

            branch.colour = colors.branchColor;
            branch.labelStyle = branch.labelStyle || {};
            branch.internalLabelStyle = branch.internalLabelStyle || {};
            branch.internalLabelStyle.colour = colors.branchLengthColor;

            if (branch.leaf) {
                branch.labelStyle.colour = colors.leafColor;
                branch.leafStyle = branch.leafStyle || {};
                branch.leafStyle.fillStyle = colors.leafColor;
                branch.leafStyle.strokeStyle = colors.leafColor;
            }
        });
    }

    newickViewer.draw();
};

window.setNewickTreeLayout = function(layout, btn) {
    if (!newickViewer) {
        return;
    }

    var targetTreeType = layout;
    var orientation = 'normal';

    if (layout === 'hierarchical-horizontal') {
        targetTreeType = 'rectangular';
    } else if (layout === 'hierarchical-vertical') {
        targetTreeType = 'hierarchical';
    } else if (layout === 'diagonal-horizontal') {
        targetTreeType = 'diagonal';
    } else if (layout === 'diagonal-vertical') {
        targetTreeType = 'diagonal';
        orientation = 'vertical';
    }

    newickViewer.setTreeType(targetTreeType);

    var viewerContainer = document.getElementById('newick-viewer');
    if (viewerContainer) {
        viewerContainer.classList.remove('orientation-diagonal-vertical');
        if (orientation === 'vertical') {
            viewerContainer.classList.add('orientation-diagonal-vertical');
        }
    }

    var layoutButtons = document.querySelectorAll('.newick-button-row .newick-mini-button');
    layoutButtons.forEach(function(button) {
        button.classList.remove('is-active');
    });

    if (btn) {
        btn.classList.add('is-active');
    }
};

window.setNewickInteraction = function(enabled) {
    var container = document.getElementById('newick-viewer');
    if (!container) {
        return;
    }

    var canvas = container.querySelector('canvas');
    if (!canvas) {
        return;
    }

    if (enabled) {
        canvas.style.pointerEvents = 'auto';
        canvas.style.cursor = 'grab';
    } else {
        canvas.style.pointerEvents = 'none';
        canvas.style.cursor = 'default';
    }
};

window.syncNewickWorkspaceHeights = function() {
    var workspace = document.querySelector('#newick_tree_viewer .newick-workspace');
    if (!workspace) {
        return;
    }

    var viewer = workspace.querySelector('#newick-viewer');
    var menu = workspace.querySelector('.newick-control-menu--side');
    if (!viewer || !menu || menu.offsetHeight === 0) {
        return;
    }

    viewer.style.height = menu.offsetHeight + 'px';

    if (newickViewer && typeof newickViewer.resizeToContainer === 'function') {
        newickViewer.resizeToContainer();
        newickViewer.draw();
    }
};

window.clearNewickForm = function() {
    var form = document.getElementById('newick-form');
    if (form) {
        form.reset();
    }

    var newickTextarea = form ? form.querySelector('textarea[name="newick"]') : null;
    if (newickTextarea) {
        newickTextarea.value = '';
    }

    newickRawTreeString = null;

    var viewerContainer = document.getElementById('newick-viewer');
    if (viewerContainer) {
        viewerContainer.innerHTML = "<div style='padding: 20px; color: #64748b;'>Paste or upload a Newick file, then click Submit to render the tree.</div>";
    }

    var interactionToggle = document.getElementById('newick-enable-interaction');
    if (interactionToggle) {
        interactionToggle.checked = false;
    }

    if (newickViewer) {
        try {
            newickViewer = null;
        } catch (err) {
            console.error('Failed to clear Newick viewer:', err);
        }
    }

    updateNewickAnnotations();
};

window.resetNewickView = function() {
    if (!newickViewer || !newickRawTreeString) {
        return;
    }

    try {
        newickViewer.setTreeType('rectangular');
        newickViewer.alignLabels = true;
        newickViewer.showBranchLengthLabels = false;
        newickViewer.showInternalNodeLabels = false;

        var viewerContainer = document.getElementById('newick-viewer');
        if (viewerContainer) {
            viewerContainer.classList.remove('orientation-diagonal-vertical');
        }

        var layoutButtons = document.querySelectorAll('.newick-button-row .newick-mini-button');
        layoutButtons.forEach(function(button) {
            button.classList.remove('is-active');
            if (button.dataset.layout === 'hierarchical-horizontal') {
                button.classList.add('is-active');
            }
        });

        if (typeof newickViewer.fitInPanel === 'function') {
            newickViewer.fitInPanel();
        } else {
            // Fallback: reload tree to reset zoom/pan state on older builds.
            newickViewer.load(newickRawTreeString);
        }

        updateNewickStyles();
        updateNewickAnnotations();
        newickViewer.draw();
    } catch (err) {
        console.error('Failed to reset Newick view:', err);
    }
};

window.downloadNewickTreeImage = function() {
    if (!newickViewer || !newickRawTreeString) {
        return;
    }

    var hiddenContainer = document.createElement('div');
    hiddenContainer.id = "hidden-newick-export-container";
    hiddenContainer.style.width = "3000px";
    hiddenContainer.style.height = "3000px";
    hiddenContainer.style.position = "absolute";
    hiddenContainer.style.left = "-9999px";
    document.body.appendChild(hiddenContainer);

    var exportViewer;
    if (Phylocanvas.default && typeof Phylocanvas.default.createTree === 'function') {
        exportViewer = Phylocanvas.default.createTree('hidden-newick-export-container');
    } else if (typeof Phylocanvas.createTree === 'function') {
        exportViewer = Phylocanvas.createTree('hidden-newick-export-container');
    } else {
        exportViewer = new Phylocanvas.Tree('hidden-newick-export-container');
    }

    exportViewer.setTreeType(newickViewer.treeType);
    exportViewer.alignLabels = newickViewer.alignLabels;
    exportViewer.showBranchLengthLabels = newickViewer.showBranchLengthLabels;
    exportViewer.showInternalNodeLabels = newickViewer.showInternalNodeLabels;

    exportViewer.textSize = 28;
    exportViewer.lineWidth = 4;

    var colors = getNewickStyleValues();
    exportViewer.customLabelColours = true;
    exportViewer.newickBranchLengthColour = colors.branchLengthColor;
    exportViewer.newickInternalNodeColour = colors.internalNodeColor;
    exportViewer.branchColour = colors.branchColor;
    exportViewer.internalLabelStyle = exportViewer.internalLabelStyle || {};
    exportViewer.internalLabelStyle.colour = colors.branchLengthColor;

    exportViewer.load(newickRawTreeString);

    if (exportViewer.branches) {
        Object.keys(exportViewer.branches).forEach(function(branchId) {
            var branch = exportViewer.branches[branchId];
            if (!branch) {
                return;
            }
            branch.colour = colors.branchColor;
            branch.labelStyle = branch.labelStyle || {};
            branch.internalLabelStyle = branch.internalLabelStyle || {};
            branch.internalLabelStyle.colour = colors.branchLengthColor;
            if (branch.leaf) {
                branch.labelStyle.colour = colors.leafColor;
                branch.leafStyle = branch.leafStyle || {};
                branch.leafStyle.fillStyle = colors.leafColor;
                branch.leafStyle.strokeStyle = colors.leafColor;
            }
        });
    }

    exportViewer.draw();

    setTimeout(function() {
        var canvas = hiddenContainer.querySelector('canvas');
        if (canvas) {
            var tempCanvas = buildAnnotatedExportCanvas(canvas, 'newick');

            var link = document.createElement('a');
            link.download = 'PhyloDendron_Newick_Tree.png';
            link.href = tempCanvas.toDataURL('image/png', 1.0);
            link.click();
        }

        document.body.removeChild(hiddenContainer);
    }, 250);
};

var inferenceViewers = {};
var inferenceRawTreeStrings = {
    bayes: PHYLODENDRON_CONFIG.bayesNewick || null,
    iqtree: PHYLODENDRON_CONFIG.iqtreeNewick || null,
    mpboot: PHYLODENDRON_CONFIG.mpbootNewick || null,
    distance: PHYLODENDRON_CONFIG.distanceNewick || null
};

function createPhylocanvasTree(containerId) {
    if (Phylocanvas.default && typeof Phylocanvas.default.createTree === 'function') {
        return Phylocanvas.default.createTree(containerId);
    }
    if (typeof Phylocanvas.createTree === 'function') {
        return Phylocanvas.createTree(containerId);
    }
    return new Phylocanvas.Tree(containerId);
}

function getInferenceStyleValues(viewerKey) {
    return {
        branchColor: (document.getElementById(viewerKey + '-branch-color') || {}).value || '#111827',
        leafColor: (document.getElementById(viewerKey + '-leaf-color') || {}).value || '#111827',
        branchLengthColor: (document.getElementById(viewerKey + '-branch-length-color') || {}).value || '#111827',
        internalNodeColor: (document.getElementById(viewerKey + '-internal-node-color') || {}).value || '#111827'
    };
}

window.updateInferenceAnnotations = function(viewerKey) {
    setTreeAnnotations(viewerKey);
};

window.updateInferenceStyles = function(viewerKey) {
    var viewer = inferenceViewers[viewerKey];
    if (!viewer) {
        return;
    }

    var colors = getInferenceStyleValues(viewerKey);
    viewer.branchColour = colors.branchColor;
    viewer.customLabelColours = true;
    viewer.newickBranchLengthColour = colors.branchLengthColor;
    viewer.newickInternalNodeColour = colors.internalNodeColor;
    viewer.internalLabelStyle = viewer.internalLabelStyle || {};
    viewer.internalLabelStyle.colour = colors.branchLengthColor;

    if (viewer.branches) {
        Object.keys(viewer.branches).forEach(function(branchId) {
            var branch = viewer.branches[branchId];
            if (!branch) {
                return;
            }

            branch.colour = colors.branchColor;
            branch.labelStyle = branch.labelStyle || {};
            branch.internalLabelStyle = branch.internalLabelStyle || {};
            branch.internalLabelStyle.colour = colors.branchLengthColor;

            if (branch.leaf) {
                branch.labelStyle.colour = colors.leafColor;
                branch.leafStyle = branch.leafStyle || {};
                branch.leafStyle.fillStyle = colors.leafColor;
                branch.leafStyle.strokeStyle = colors.leafColor;
            }
        });
    }

    viewer.draw();
};

window.setInferenceTreeLayout = function(viewerKey, layout, btn) {
    var viewer = inferenceViewers[viewerKey];
    if (!viewer) {
        return;
    }

    var targetTreeType = layout;
    var orientation = 'normal';

    if (layout === 'hierarchical-horizontal') {
        targetTreeType = 'rectangular';
    } else if (layout === 'hierarchical-vertical') {
        targetTreeType = 'hierarchical';
    } else if (layout === 'diagonal-horizontal') {
        targetTreeType = 'diagonal';
    } else if (layout === 'diagonal-vertical') {
        targetTreeType = 'diagonal';
        orientation = 'vertical';
    }

    viewer.setTreeType(targetTreeType);

    var viewerContainer = document.getElementById(viewerKey + '-viewer');
    if (viewerContainer) {
        viewerContainer.classList.remove('orientation-diagonal-vertical');
        if (orientation === 'vertical') {
            viewerContainer.classList.add('orientation-diagonal-vertical');
        }
    }

    var workspace = document.getElementById(viewerKey + '-workspace');
    if (workspace) {
        var layoutButtons = workspace.querySelectorAll('.newick-button-row .newick-mini-button');
        layoutButtons.forEach(function(button) {
            button.classList.remove('is-active');
        });
        if (btn) {
            btn.classList.add('is-active');
        }
    }

    viewer.draw();
};

window.setInferenceAlignLabels = function(viewerKey, enabled) {
    var viewer = inferenceViewers[viewerKey];
    if (!viewer) {
        return;
    }
    viewer.alignLabels = enabled;
    viewer.draw();
};

window.setInferenceBranchLengths = function(viewerKey, enabled) {
    var viewer = inferenceViewers[viewerKey];
    if (!viewer) {
        return;
    }
    viewer.showBranchLengthLabels = enabled;
    updateInferenceStyles(viewerKey);
};

window.setInferenceInternalNodes = function(viewerKey, enabled) {
    var viewer = inferenceViewers[viewerKey];
    if (!viewer) {
        return;
    }
    viewer.showInternalNodeLabels = enabled;
    updateInferenceStyles(viewerKey);
};

window.setInferenceInteraction = function(viewerKey, enabled) {
    var container = document.getElementById(viewerKey + '-viewer');
    if (!container) {
        return;
    }

    var canvas = container.querySelector('canvas');
    if (!canvas) {
        return;
    }

    if (enabled) {
        canvas.style.pointerEvents = 'auto';
        canvas.style.cursor = 'grab';
    } else {
        canvas.style.pointerEvents = 'none';
        canvas.style.cursor = 'default';
    }
};

window.syncInferenceWorkspaceHeights = function(viewerKey) {
    var workspace = document.getElementById(viewerKey + '-workspace');
    if (!workspace) {
        return;
    }

    var viewer = document.getElementById(viewerKey + '-viewer');
    var menu = workspace.querySelector('.newick-control-menu--side');
    if (!viewer || !menu || menu.offsetHeight === 0) {
        return;
    }

    viewer.style.height = menu.offsetHeight + 'px';

    var treeViewer = inferenceViewers[viewerKey];
    if (treeViewer && typeof treeViewer.resizeToContainer === 'function') {
        treeViewer.resizeToContainer();
        treeViewer.draw();
    }
};

window.resetInferenceView = function(viewerKey) {
    var viewer = inferenceViewers[viewerKey];
    var rawTree = inferenceRawTreeStrings[viewerKey];
    if (!viewer || !rawTree) {
        return;
    }

    try {
        viewer.setTreeType('rectangular');
        viewer.alignLabels = true;
        viewer.showBranchLengthLabels = false;
        viewer.showInternalNodeLabels = false;

        var viewerContainer = document.getElementById(viewerKey + '-viewer');
        if (viewerContainer) {
            viewerContainer.classList.remove('orientation-diagonal-vertical');
        }

        var workspace = document.getElementById(viewerKey + '-workspace');
        if (workspace) {
            var layoutButtons = workspace.querySelectorAll('.newick-button-row .newick-mini-button');
            layoutButtons.forEach(function(button) {
                button.classList.remove('is-active');
                if (button.dataset.layout === 'hierarchical-horizontal') {
                    button.classList.add('is-active');
                }
            });
        }

        if (typeof viewer.fitInPanel === 'function') {
            viewer.fitInPanel();
        } else {
            viewer.load(rawTree);
        }

        updateInferenceStyles(viewerKey);
        updateInferenceAnnotations(viewerKey);

        var interactionToggle = document.getElementById(viewerKey + '-enable-interaction');
        if (interactionToggle) {
            interactionToggle.checked = false;
        }
        setInferenceInteraction(viewerKey, false);
    } catch (err) {
        console.error('Failed to reset inference view for ' + viewerKey + ':', err);
    }
};

window.downloadInferenceTreeImage = function(viewerKey, fileName) {
    var viewer = inferenceViewers[viewerKey];
    var rawTree = inferenceRawTreeStrings[viewerKey];
    if (!viewer || !rawTree) {
        return;
    }

    var hiddenContainerId = 'hidden-' + viewerKey + '-export-container';
    var hiddenContainer = document.createElement('div');
    hiddenContainer.id = hiddenContainerId;
    hiddenContainer.style.width = '3000px';
    hiddenContainer.style.height = '3000px';
    hiddenContainer.style.position = 'absolute';
    hiddenContainer.style.left = '-9999px';
    document.body.appendChild(hiddenContainer);

    var exportViewer = createPhylocanvasTree(hiddenContainerId);
    exportViewer.setTreeType(viewer.treeType);
    exportViewer.alignLabels = viewer.alignLabels;
    exportViewer.showBranchLengthLabels = viewer.showBranchLengthLabels;
    exportViewer.showInternalNodeLabels = viewer.showInternalNodeLabels;
    exportViewer.textSize = 28;
    exportViewer.lineWidth = 4;

    var colors = getInferenceStyleValues(viewerKey);
    exportViewer.customLabelColours = true;
    exportViewer.newickBranchLengthColour = colors.branchLengthColor;
    exportViewer.newickInternalNodeColour = colors.internalNodeColor;
    exportViewer.branchColour = colors.branchColor;
    exportViewer.internalLabelStyle = exportViewer.internalLabelStyle || {};
    exportViewer.internalLabelStyle.colour = colors.branchLengthColor;
    exportViewer.load(rawTree);

    if (exportViewer.branches) {
        Object.keys(exportViewer.branches).forEach(function(branchId) {
            var branch = exportViewer.branches[branchId];
            if (!branch) {
                return;
            }
            branch.colour = colors.branchColor;
            branch.labelStyle = branch.labelStyle || {};
            branch.internalLabelStyle = branch.internalLabelStyle || {};
            branch.internalLabelStyle.colour = colors.branchLengthColor;
            if (branch.leaf) {
                branch.labelStyle.colour = colors.leafColor;
                branch.leafStyle = branch.leafStyle || {};
                branch.leafStyle.fillStyle = colors.leafColor;
                branch.leafStyle.strokeStyle = colors.leafColor;
            }
        });
    }

    exportViewer.draw();

    setTimeout(function() {
        var canvas = hiddenContainer.querySelector('canvas');
        if (canvas) {
            var tempCanvas = buildAnnotatedExportCanvas(canvas, viewerKey);

            var link = document.createElement('a');
            link.download = fileName || ('PhyloDendron_' + viewerKey + '_Tree.png');
            link.href = tempCanvas.toDataURL('image/png', 1.0);
            link.click();
        }

        document.body.removeChild(hiddenContainer);
    }, 250);
};

function initInferenceViewer(viewerKey) {
    var rawTree = inferenceRawTreeStrings[viewerKey];
    if (!rawTree) {
        return;
    }

    var viewerContainerId = viewerKey + '-viewer';
    function tryInitialize() {
        var viewerContainer = document.getElementById(viewerContainerId);
        if (!viewerContainer) return;

        // wait until the container is visible
        if (viewerContainer.offsetWidth === 0) {
            setTimeout(tryInitialize, 100);
            return;
        }

        try {
            // prevent initializing twice
            if (!inferenceViewers[viewerKey]) {
                var viewer = createPhylocanvasTree(viewerContainerId);
                viewer.setTreeType('rectangular');
                viewer.alignLabels = true;
                viewer.textSize = 14;
                viewer.lineWidth = 2;
                viewer.load(rawTree);
                inferenceViewers[viewerKey] = viewer;

                updateInferenceStyles(viewerKey);
                updateInferenceAnnotations(viewerKey);
                syncInferenceWorkspaceHeights(viewerKey);
                setInferenceInteraction(viewerKey, false);

                var interactionToggle = document.getElementById(viewerKey + '-enable-interaction');
                if (interactionToggle) {
                    interactionToggle.checked = false;
                }
            }
        } catch (err) {
            viewerContainer.innerHTML = "<div style='padding: 20px; color: red; font-weight: bold;'>Error rendering tree: " + err.message + "</div>";
        }
    }
    
    // start the process
    tryInitialize();
}

window.onload = function() {
    if (PHYLODENDRON_CONFIG.activeTab) {
        document.querySelectorAll('.tool').forEach(function(t) {
            t.classList.remove('active');
        });
        var activeTool = document.getElementById(PHYLODENDRON_CONFIG.activeTab);
        if (activeTool) {
            activeTool.classList.add('active');
        }
    }

    if (PHYLODENDRON_CONFIG.autoDownloadUrl) {
        setTimeout(function() {
            var autoLink = document.createElement('a');
            autoLink.href = PHYLODENDRON_CONFIG.autoDownloadUrl;
            autoLink.setAttribute('download', '');
            document.body.appendChild(autoLink);
            autoLink.click();
            autoLink.remove();
        }, 150);
    }

    if (newickRawTreeString) {
        function initNewick() {
            var newickContainer = document.getElementById('newick-viewer');
            if (!newickContainer) return;

            // wait until the container is visible
            if (newickContainer.offsetWidth === 0) {
                setTimeout(initNewick, 100);
                return;
            }

            try {
                enableSplitNewickLabelColors();

                if (Phylocanvas.default && typeof Phylocanvas.default.createTree === 'function') {
                    newickViewer = Phylocanvas.default.createTree('newick-viewer');
                } else if (typeof Phylocanvas.createTree === 'function') {
                    newickViewer = Phylocanvas.createTree('newick-viewer');
                } else {
                    newickViewer = new Phylocanvas.Tree('newick-viewer');
                }

                newickViewer.setTreeType('rectangular');
                newickViewer.alignLabels = true;
                newickViewer.textSize = 14;
                newickViewer.lineWidth = 2;

                newickViewer.load(newickRawTreeString);
                updateNewickStyles();
                updateNewickAnnotations();

                syncNewickWorkspaceHeights();
                window.addEventListener('resize', syncNewickWorkspaceHeights);

                // Default to scroll-friendly mode
                setNewickInteraction(false);
                var newickInteractionToggle = document.getElementById('newick-enable-interaction');
                if (newickInteractionToggle) {
                    newickInteractionToggle.checked = false;
                }
            } catch (err) {
                newickContainer.innerHTML = "<div style='padding: 20px; color: red; font-weight: bold;'>Error rendering Newick tree: " + err.message + "</div>";
            }
        }
        
        // start the process
        initNewick();
    }

    ['bayes', 'iqtree', 'mpboot', 'distance'].forEach(function(viewerKey) {
        initInferenceViewer(viewerKey);
    });

    window.addEventListener('resize', function() {
        ['bayes', 'iqtree', 'mpboot', 'distance'].forEach(function(viewerKey) {
            syncInferenceWorkspaceHeights(viewerKey);
        });
    });

    // Require molecule type only for formats that need it.
    var conversionForm = document.getElementById("conversion-form");
    var conversionOutputFormat = document.getElementById("conversion-output-format");
    var moleculeRadios = document.querySelectorAll(".conversion-molecule-type");
    var moleculeHelp = document.getElementById("conversion-molecule-help");

    moleculeRadios.forEach(function(radio) {
        radio.dataset.wasChecked = radio.checked ? "true" : "false";
        radio.addEventListener("click", function() {
            if (this.dataset.wasChecked === "true") {
                this.checked = false;
            }

            moleculeRadios.forEach(function(otherRadio) {
                otherRadio.dataset.wasChecked = otherRadio.checked ? "true" : "false";
            });
        });
    });

    function updateConversionMoleculeRequirement() {
        if (!conversionOutputFormat || moleculeRadios.length === 0) {
            return;
        }

        var requiresMoleculeType = ["nexus", "genbank", "embl"].indexOf(conversionOutputFormat.value) !== -1;

        moleculeRadios.forEach(function(radio) {
            radio.required = requiresMoleculeType;
            radio.disabled = !requiresMoleculeType;
            radio.setCustomValidity("");

            if (!requiresMoleculeType) {
                radio.checked = false;
            }
        });

        if (!requiresMoleculeType) {
            moleculeRadios.forEach(function(radio) {
                radio.dataset.wasChecked = "false";
            });
        }

        if (moleculeHelp) {
            moleculeHelp.textContent = requiresMoleculeType
                ? "Required for selected output format."
                : "";
        }
    }

    if (conversionOutputFormat) {
        conversionOutputFormat.addEventListener("change", updateConversionMoleculeRequirement);
    }

    if (conversionForm) {
        conversionForm.addEventListener("submit", function() {
            updateConversionMoleculeRequirement();
        });
    }

    updateConversionMoleculeRequirement();

    // ----------------------------------------------------
    // LOADING SPINNER LOGIC (Grouped for all tools)
    // ----------------------------------------------------
    var asyncForms = [
        "#distance-form", 
        "#mafft-form", 
        "#clustalo-form", 
        "#muscle-form",
        "#parsimony-form",
        "#iqtree-form"
    ];

    asyncForms.forEach(function(formId) {
        var form = document.querySelector(formId);
        if (form) {
            form.addEventListener("submit", function() {
                document.getElementById("loading-overlay").style.display = "flex";
            });
        }
    });

    // The Bayesian endpoint returns a downloadable file; use fetch so we can hide the overlay
    // only after the browser has received the response.
    var bayesForm = document.getElementById("bayesian-inference-form");
    if (bayesForm) {
        bayesForm.addEventListener("submit", async function(event) {
            event.preventDefault();

            var overlay = document.getElementById("loading-overlay");
            overlay.style.display = "flex";

            try {
                var formData = new FormData(bayesForm);
                var response = await fetch(bayesForm.action, {
                    method: "POST",
                    body: formData
                });

                if (!response.ok) {
                    throw new Error("Request failed with status " + response.status);
                }

                var contentType = (response.headers.get("Content-Type") || "").toLowerCase();
                if (contentType.indexOf("text/html") !== -1) {
                    // Render server-side validation/errors returned as HTML.
                    var html = await response.text();
                    document.open();
                    document.write(html);
                    document.close();
                    return;
                }

                var blob = await response.blob();
                var disposition = response.headers.get("Content-Disposition") || "";
                var fallbackName = "mrbayes_results.zip";
                var fileName = fallbackName;

                var utf8NameMatch = disposition.match(/filename\*=UTF-8''([^;]+)/i);
                var asciiNameMatch = disposition.match(/filename=\"?([^\";]+)\"?/i);
                if (utf8NameMatch && utf8NameMatch[1]) {
                    fileName = decodeURIComponent(utf8NameMatch[1]);
                } else if (asciiNameMatch && asciiNameMatch[1]) {
                    fileName = asciiNameMatch[1];
                }

                var downloadUrl = URL.createObjectURL(blob);
                var link = document.createElement("a");
                link.href = downloadUrl;
                link.download = fileName;
                document.body.appendChild(link);
                link.click();
                link.remove();
                URL.revokeObjectURL(downloadUrl);
            } catch (error) {
                alert("Bayesian inference failed: " + error.message);
            } finally {
                overlay.style.display = "none";
            }
        });
    }

    // ----------------------------------------------------
    // MRBAYES PRESET LOGIC
    // ----------------------------------------------------
    window.applyMrBayesPreset = function() {
        var preset = document.getElementById("mrbayes-preset").value;
        var ngen = document.getElementById("mrbayes-ngen");
        var samplefreq = document.getElementById("mrbayes-samplefreq");
        var printfreq = document.getElementById("mrbayes-printfreq");
        var burnin = document.getElementById("mrbayes-burnin");

        if (preset === "intensif") {
            ngen.value = 1000000;
            samplefreq.value = 100;
            printfreq.value = 100000;
            burnin.value = 2500;
        } else if (preset === "rapide") {
            ngen.value = 100000;
            samplefreq.value = 100;
            printfreq.value = 10000;
            burnin.value = 250;
        }
    };

    // MrBayes Custom Dropdown Trigger
    ["mrbayes-ngen", "mrbayes-samplefreq", "mrbayes-printfreq", "mrbayes-burnin"].forEach(function(id) {
        var el = document.getElementById(id);
        if(el) el.addEventListener('input', () => document.getElementById("mrbayes-preset").value = "custom");
    });

    // MPBOOT & IQTREE PRESET LOGIC
    // ----------------------------------------------------
    
    // MPBoot Presets
    window.applyMPBootPreset = function() {
        var preset = document.getElementById("mpboot-preset").value;
        var bb = document.getElementById("mpb-bb");

        if (preset === "standard") {
            bb.value = 1000; 
        } else if (preset === "fast") {
            bb.value = 0; // Skips bootstraps for a quick tree
        } else if (preset === "thorough") {
            bb.value = 2000; 
        }
    };

    // MPBoot Custom Dropdown Trigger
    var mpbBb = document.getElementById("mpb-bb");
    if(mpbBb) {
        mpbBb.addEventListener('input', () => document.getElementById("mpboot-preset").value = "custom");
    };

    // IQ-TREE Presets
    window.applyIQTreePreset = function() {
        var preset = document.getElementById("iqtree-preset").value;
        var m = document.getElementById("iq-m");
        var bb = document.getElementById("iq-bb");
        var alrt = document.getElementById("iq-alrt");

        if (preset === "standard") {
            m.value = "MFP"; bb.value = 1000; alrt.value = 1000;
        } else if (preset === "fast") {
            m.value = "MFP"; bb.value = 0; alrt.value = 0;
        } else if (preset === "thorough") {
            m.value = "MFP"; bb.value = 2000; alrt.value = 2000; 
        }
    };

    // IQ-TREE Custom Dropdown Trigger
    ["iq-m", "iq-bb", "iq-alrt"].forEach(function(id) {
        var el = document.getElementById(id);
        if(el) el.addEventListener('input', () => document.getElementById("iqtree-preset").value = "custom");
    });
};
