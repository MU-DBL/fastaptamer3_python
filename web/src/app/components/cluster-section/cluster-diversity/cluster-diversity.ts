import { Component, signal } from '@angular/core';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';

@Component({
  selector: 'app-cluster-diversity',
  imports: [   
    CommonModule,
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster-diversity.html',
  styleUrl: './cluster-diversity.scss',
})
export class ClusterDiversity {

  adjustClusterMetadata: string = 'no';
  clusterMetaXAxis:  string = "Cluster";
  clusterMetaYAxis1:  string ="Seq. count";
  clusterMetaColor1: string = "#0000ff";
  clusterMetaYAxis2:  string ="Read count";
  clusterMetaColor2:  string ="#ffa500";
  clusterMetaYAxis3:  string ="Avg. LED";
  clusterMetaColor3:  string ="#228b22";
  clusterMetaPlotTitle:  string ="Cluster metaplots";

  
  kmerSize:  string ="3";
  plotType:  string ="UMAP";

  adjustKmerPlot: string = 'no';

  kmerPlotXAxis:  string ="Dim1";
  kmerPlotYAxis:  string ="Dim2";
  kmerPlotLegendTitle:  string ="Cluster";
  kmerPlotTitle:  string ="Cluster k-mer plot";
  kmerPlotColorPalette:  string ="#000000";

  clusterMetaPlot() {
      throw new Error('Method not implemented.');
  }

  kmerPlot() {
    throw new Error('Method not implemented.');
  }
  

    onDownload() {
    throw new Error('Method not implemented.');
    }
    onStart() {
    throw new Error('Method not implemented.');
    }

    isProcessing = signal(false);
    processedFileName = signal('');
    
    selectedFile: File | null = null;
    savedFileName: string = '';
    uploadComplete: boolean = false;

    number_reads_to_cluster:number = 70
    led:number = 19
    number_of_cluster:number = 20


    onFileSelected(result: FileUploadResult): void {
      this.selectedFile = result.file;
      this.processedFileName.set('');
      console.log('File 1 selected:', result.fileName);
    }
  
    onUploadComplete(result: FileUploadResult): void {
      if (result.uploadComplete && result.savedFileName) {
        this.uploadComplete = true;
        this.savedFileName = result.savedFileName;
        console.log('Upload 1 complete:', result.savedFileName);
      } else if (result.error) {
        console.error('Upload 1 failed:', result.error);
      }
    }
}
