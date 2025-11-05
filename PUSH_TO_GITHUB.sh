#!/bin/bash
# Script to push LassaPartitioner to GitHub

echo "=== LassaPartitioner GitHub Push Script ==="
echo ""
echo "This script will help you push LassaPartitioner to GitHub."
echo "Make sure you have:"
echo "  1. Created a new repository on GitHub named 'lassa_partitioner'"
echo "  2. Your GitHub username ready"
echo ""
read -p "Enter your GitHub username: " GITHUB_USERNAME

echo ""
echo "Adding all files..."
git add .

echo "Making initial commit..."
git commit -m "Initial commit: LassaPartitioner v0.1.0 - Tool for creating APOBEC3 partitions from Lassa virus alignments"

echo "Renaming branch to main..."
git branch -M main

echo "Adding remote origin..."
git remote add origin https://github.com/${GITHUB_USERNAME}/lassa_partitioner.git

echo "Pushing to GitHub..."
git push -u origin main

echo ""
echo "Done! Your repository should now be available at:"
echo "https://github.com/${GITHUB_USERNAME}/lassa_partitioner"
