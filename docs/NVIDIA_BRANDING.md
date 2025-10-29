# NVIDIA Branding Implementation Guide

## 🎨 NVIDIA Brand Colors

### Primary Colors
- **NVIDIA Green**: `#76B900`
  - Primary brand color
  - Used for: Headers, primary buttons, success states, highlights
  
- **NVIDIA Dark**: `#1A1A1A`
  - Secondary brand color
  - Used for: Text, sidebar backgrounds, headers

### Secondary Colors
- **NVIDIA Accent (Teal)**: `#00D4AA`
  - Used for: Secondary buttons, info messages, accents
  
- **NVIDIA Gray**: `#666666`
  - Used for: Body text, secondary information
  
- **NVIDIA Light Gray**: `#E5E5E5`
  - Used for: Borders, backgrounds, dividers

### Status Colors
- **Success/Complete**: `#76B900` (NVIDIA Green)
- **Active/In-Progress**: `#00D4AA` (Accent Teal)
- **Warning**: `#FFA500` (Orange)
- **Error**: `#FF3838` (Red)
- **Pending**: `#666666` (Gray)

## 🖼️ Logo Usage

### Logo Location
- **Path**: `image/nvidia.jpg`
- **Placement**: 
  - Top-left of main content (120px width)
  - Sidebar header (full width)

### Implementation
```python
# Main header logo
st.image("image/nvidia.jpg", width=120)

# Sidebar logo
st.sidebar.image("image/nvidia.jpg", use_container_width=True)
```

## 🎯 Design Patterns

### Headers
```python
st.markdown("""
<div style="background: linear-gradient(135deg, #1A1A1A 0%, #2D2D2D 100%); 
            padding: 20px; border-radius: 8px; border-left: 6px solid #76B900;">
    <h1 style="color: #FFFFFF; margin: 0; border: none;">
        🧬 Title Here
    </h1>
    <p style="color: #76B900; font-weight: 600; margin: 10px 0 0 0;">
        Subtitle or tagline
    </p>
</div>
""", unsafe_allow_html=True)
```

### Buttons
- **Primary**: Green gradient (`#76B900` → `#5A8F00`)
- **Secondary**: Gray gradient (`#666666` → `#4D4D4D`)
- **Download**: Teal gradient (`#00D4AA` → `#00A88A`)

### Progress Indicators
- **Progress Bar**: Green-to-Teal gradient
- **Stage Complete**: Green with glow effect
- **Stage Active**: Teal with pulse animation
- **Stage Pending**: Gray

### Cards & Containers
```python
# Info Card
st.markdown("""
<div style="background-color: #f8f9fa; padding: 15px; 
            border-radius: 6px; border-left: 4px solid #00D4AA;">
    <p>Content here</p>
</div>
""", unsafe_allow_html=True)

# Metric Card
background: linear-gradient(135deg, #f8f9fa 0%, #e5e5e5 100%);
border-left: 4px solid #76B900;
box-shadow: 0 2px 8px rgba(0,0,0,0.1);
```

## 📱 Application Layouts

### Workflow App (binding_workflow_app.py)

#### Header Layout
```
┌─────────────────────────────────────────┐
│  [LOGO]  NVIDIA Protein Binding Design  │
│          ESMFold → RFDiffusion → ...     │
└─────────────────────────────────────────┘
```

#### Sidebar
```
┌──────────────────┐
│   [NVIDIA LOGO]  │
│                  │
│ 🧬 Binding Wflow │
│ ─────────────── │
│ 📁 Project Info  │
│ ⚙️ Configuration │
│ 💾 Session Mgmt  │
└──────────────────┘
```

### Single Structure App (app_v2.py)

#### Header Layout
```
┌─────────────────────────────────────────┐
│  [LOGO]  NVIDIA Protein Structure Pred  │
│          Powered by NVIDIA AI Models     │
└─────────────────────────────────────────┘
```

## 🎨 CSS Classes

### Stage Status
```css
.stage-complete   /* ✅ Green with glow */
.stage-active     /* 🔄 Teal with pulse */
.stage-pending    /* ⭕ Gray */
.stage-failed     /* ❌ Red */
```

### Quality Badges
```css
.quality-excellent  /* A: Green gradient */
.quality-good      /* B: Teal gradient */
.quality-moderate  /* C: Orange gradient */
.quality-poor      /* D-F: Red gradient */
```

### NVIDIA Header
```css
.nvidia-header {
    background: linear-gradient(135deg, #1A1A1A 0%, #2D2D2D 100%);
    padding: 20px;
    border-radius: 8px;
    border-left: 6px solid #76B900;
}
```

## 📋 Typography

### Font Family
```css
font-family: 'NVIDIA Sans', 'Helvetica Neue', Arial, sans-serif;
```

### Font Weights
- **Bold Headers**: 700
- **Medium**: 600
- **Regular**: 500
- **Light**: 400

### Text Colors
- **Primary Text**: `#1A1A1A` (NVIDIA Dark)
- **Secondary Text**: `#666666` (NVIDIA Gray)
- **Light Text**: `#999999`
- **Sidebar Text**: `#E5E5E5` (on dark background)

## 🌟 Visual Effects

### Gradients
```css
/* Primary Button */
background: linear-gradient(135deg, #76B900 0%, #5A8F00 100%);

/* Progress Bar */
background: linear-gradient(90deg, #76B900 0%, #00D4AA 100%);

/* Sidebar */
background: linear-gradient(180deg, #1A1A1A 0%, #2D2D2D 100%);

/* Card */
background: linear-gradient(135deg, #f8f9fa 0%, #e5e5e5 100%);
```

### Shadows
```css
/* Card Shadow */
box-shadow: 0 2px 8px rgba(0,0,0,0.1);

/* Hover Shadow (Green) */
box-shadow: 0 4px 12px rgba(118, 185, 0, 0.4);

/* Button Shadow */
box-shadow: 0 2px 8px rgba(118, 185, 0, 0.3);
```

### Animations
```css
/* Pulse for active stage */
@keyframes pulse {
    0%, 100% { opacity: 1; }
    50% { opacity: 0.6; }
}

/* Hover lift */
transition: transform 0.2s;
transform: translateY(-2px);
```

## 📐 Spacing & Sizing

### Border Radius
- **Large cards**: `8px`
- **Buttons**: `6px`
- **Small elements**: `4px`

### Padding
- **Headers**: `20px`
- **Cards**: `15-20px`
- **Buttons**: `10px 24px`

### Border Width
- **Accent left border**: `4-6px`
- **Regular borders**: `1-2px`

## 🎯 Implementation Checklist

### Both Apps
- ✅ NVIDIA logo in header
- ✅ NVIDIA logo in sidebar
- ✅ Green primary color scheme
- ✅ Dark sidebar background
- ✅ Branded footer with copyright
- ✅ Gradient buttons
- ✅ Success messages with green
- ✅ Proper spacing and shadows

### Workflow App Specific
- ✅ Stage progress with color coding
- ✅ Quality score badges
- ✅ Metric cards with green accent
- ✅ Pipeline visualization

### Single Structure App Specific
- ✅ Model info cards
- ✅ Confidence score highlighting
- ✅ PDB viewer container

## 📊 Contrast Ratios (Accessibility)

- **NVIDIA Green (#76B900) on White**: ✅ AA compliant
- **NVIDIA Dark (#1A1A1A) on White**: ✅ AAA compliant
- **White on NVIDIA Dark**: ✅ AAA compliant
- **NVIDIA Green on Dark**: ✅ AA compliant

## 🔗 Brand Resources

### Official NVIDIA Brand Guidelines
- NVIDIA brand colors are official
- Maintain 6px minimum border for logo
- Never distort or modify logo
- Use official green (#76B900) consistently

### Application-Specific
- **Workflow App**: Green = success/complete
- **Single App**: Green = high confidence
- **Both**: Teal = secondary actions/info
- **Both**: Dark sidebar for controls

## 💡 Usage Examples

### Success Message
```python
st.success("✅ Structure prediction completed!")
# Renders with green background and green border
```

### Primary Button
```python
st.button("Predict Structure", type="primary")
# Renders with green gradient background
```

### Info Card
```python
st.info("ℹ️ Tip: Use OpenFold2 for faster predictions")
# Renders with teal border
```

### Metric Display
```python
st.metric("Quality Score", "85/100")
# Value renders in NVIDIA green
```

## 🚀 Future Enhancements

### Planned
- [ ] Custom NVIDIA font import
- [ ] Loading spinner with NVIDIA branding
- [ ] Progress animations
- [ ] Interactive logo on hover
- [ ] Dark mode toggle (maintain green accent)
- [ ] Accessibility improvements
- [ ] Mobile-responsive adjustments

---

**Note**: All styling follows NVIDIA's official brand guidelines while maintaining application usability and accessibility standards.
