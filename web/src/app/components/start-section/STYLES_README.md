# Shared Styles for Startpage Components

This folder contains shared styles that can be reused across all startpage components (preprocess, count, recount, etc.).

## Usage

Simply import the shared styles file at the top of your component's SCSS file:

```scss
@import '../shared-styles.scss';
```

## Available Styles

### Layout & Container
- `.container` - Main container with responsive width and padding
- `.form-section` - Section wrapper for form elements
- `.full-width` - Full width utility class

### Typography & Headings
- `.section-title` - Section headings
- `.section-label` - Form labels
- `.link` - Styled links with hover effects

### Form Elements
- `.radio-group` - Radio button group layout
- Material form field overrides for consistent styling

### Slider Components
- `.slider-group` - Slider wrapper
- `.slider-header` - Slider value display header
- `.value-chip` - Chip-style value indicators
- `.dual-range-container` - Dual range slider container
- `.range-slider` - Slider styling
- `.info-icon` - Help icon with tooltip

### Buttons
- `.action-buttons` - Button group container
- `.start-button` - Primary action button (dark gray)
- `.download-button` - Download button (dark gray, turns blue when ready)
- `.download-ready` - Applied to download button when file is ready
- `.browse-button` - File browse button

### Upload Elements
- `.file-upload-section` - Upload section layout
- `.file-name` - File name display
- `.upload-note` - Upload instruction text

### Progress Indicators
- `.progress-container` - Progress bar wrapper
- `.progress-text` - Progress status text
- `.spinning` - Spinning animation class (for loading states)

### Animations
- `spin` - Rotation animation (360° continuous)
- `pulse` - Pulsing glow animation (for download-ready state)

## Component-Specific Styles

Each component can still have its own specific styles. Just add them below the import statement:

```scss
@import '../shared-styles.scss';

// Component-specific styles
.my-custom-class {
  // your custom styles here
}
```

## Responsive Design

All shared styles include responsive breakpoints for mobile devices (max-width: 600px):
- Stacked layouts on mobile
- Full-width buttons
- Adjusted spacing and sizing

## Benefits

✅ **Consistency** - All components look and behave the same way
✅ **Maintainability** - Update styles in one place
✅ **Reusability** - Easy to add new components with consistent styling
✅ **Smaller Bundle** - No duplicate CSS in final build
✅ **Easier Updates** - Change design system-wide with minimal effort

## Example

```html
<div class="container">
  <div class="form-section">
    <h3 class="section-title">Upload File</h3>
    <!-- your content -->
  </div>
  
  <div class="action-buttons">
    <button class="start-button" [class.spinning]="isProcessing">
      Start
    </button>
    <button class="download-button" [class.download-ready]="hasFile">
      Download
    </button>
  </div>
</div>
```

## Updates

When adding new shared styles:
1. Add them to `shared-styles.scss`
2. Update this README with documentation
3. Test across all components using the shared styles
