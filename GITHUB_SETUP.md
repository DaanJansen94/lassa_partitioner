# Setting Up LassaPartitioner on GitHub

This guide will help you push LassaPartitioner to GitHub as a new repository.

## Step 1: Create a New Repository on GitHub

1. Go to https://github.com
2. Click the "+" icon in the top right corner
3. Select "New repository"
4. Repository name: `lassa_partitioner`
5. Description: "A bioinformatics tool for creating APOBEC3 and non-APOBEC3 partitions from Lassa virus sequence alignments"
6. Choose Public or Private
7. **DO NOT** initialize with README, .gitignore, or license (we already have these)
8. Click "Create repository"

## Step 2: Initialize Git in Your Local Directory

```bash
cd /data/Daan/Projects/Lassa/Lassa_mouse/ADAR/lassa_partitioner
git init
```

## Step 3: Add All Files

```bash
git add .
```

## Step 4: Make Initial Commit

```bash
git commit -m "Initial commit: LassaPartitioner v0.1.0"
```

## Step 5: Add GitHub Remote

Replace `YOUR_USERNAME` with your actual GitHub username:

```bash
git remote add origin https://github.com/YOUR_USERNAME/lassa_partitioner.git
```

Or if you prefer SSH:

```bash
git remote add origin git@github.com:YOUR_USERNAME/lassa_partitioner.git
```

## Step 6: Push to GitHub

```bash
git branch -M main
git push -u origin main
```

## Step 7: Verify

Visit `https://github.com/YOUR_USERNAME/lassa_partitioner` to verify your repository is uploaded.

## Optional: Add GitHub Pages Documentation

If you want to create a GitHub Pages site:

1. Go to Settings → Pages
2. Select "main" branch and "/docs" folder (or root)
3. Your README.md will be automatically rendered

## File Structure

Your repository should have:
```
lassa_partitioner/
├── .gitignore
├── LICENSE
├── README.md
├── setup.py
└── lassa_partitioner/
    ├── __init__.py
    ├── cli.py
    └── core.py
```

## Future Updates

To push updates in the future:

```bash
git add .
git commit -m "Description of changes"
git push origin main
```

