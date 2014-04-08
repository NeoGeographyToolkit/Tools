// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__
#include <cassert>
#include <iostream>
#include <QDialogButtonBox>
#include <QLabel>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QWidget>
#include <QHeaderView>
#include "chooseFilesDlg.h"
#include "utils.h"
using namespace std;
using namespace utils;

// Allow the user to choose which files to hide/show in the GUI.
// User's choice will be processed by aspView::showFilesChosenByUser().
chooseFilesDlg::chooseFilesDlg(QWidget * parent,
                               const std::vector<polyOptions> & optionsVec):
  QWidget(parent){

  setWindowModality(Qt::ApplicationModal); 

  int spacing = 6;
  
  QVBoxLayout * vBoxLayout = new QVBoxLayout(this);
  vBoxLayout->setSpacing(spacing);
  vBoxLayout->setAlignment(Qt::AlignLeft);

  // Label
  //QLabel * label = new QLabel(chooseFilesDlg::selectFilesTag(), this);

  // The layout having the file names. It will be filled in dynamically later.
  m_filesTable = new QTableWidget();

//   QDialogButtonBox * submitBox = new QDialogButtonBox(this);
//   submitBox->setOrientation(Qt::Horizontal);
//   submitBox->setStandardButtons(QDialogButtonBox::Ok);

  m_filesTable->horizontalHeader()->hide();
  m_filesTable->verticalHeader()->hide();
    
  //vBoxLayout->addWidget(label);
  vBoxLayout->addWidget(m_filesTable);
//   vBoxLayout->addWidget(submitBox);

//   connect(submitBox, SIGNAL(accepted()), this, SLOT(accept()));
//   connect(submitBox, SIGNAL(rejected()), this, SLOT(reject()));

  return;
}

chooseFilesDlg::~chooseFilesDlg(){}

void chooseFilesDlg::chooseFiles(const std::vector<polyOptions> & optionsVec){

  // See the top of this file for documentation.
  
  int numFiles = optionsVec.size();
  int numCols = 2;
  m_filesTable->setRowCount(numFiles);
  m_filesTable->setColumnCount(numCols);

  for (int fileIter = 0; fileIter < numFiles; fileIter++){

    QTableWidgetItem *item;    
    item = new QTableWidgetItem(1);
    item->data(Qt::CheckStateRole);
    item->setCheckState(Qt::Checked);
    m_filesTable->setItem(fileIter, 0, item);
    
    string fileName = optionsVec[fileIter].polyFileName;
    item = new QTableWidgetItem(fileName.c_str());
    item->setFlags(Qt::NoItemFlags);
    item->setForeground(QColor::fromRgb(0,0,0));
    m_filesTable->setItem(fileIter, numCols - 1, item);
    
  }
  
  QStringList rowNamesList;
  for (int fileIter = 0; fileIter < numFiles; fileIter++) rowNamesList << "";
  m_filesTable->setVerticalHeaderLabels(rowNamesList);

  QStringList colNamesList;
  for (int colIter = 0; colIter < numCols; colIter++) colNamesList << "";
  m_filesTable->setHorizontalHeaderLabels(colNamesList);
  QTableWidgetItem * hs = m_filesTable->horizontalHeaderItem(0);
  hs->setBackground(QBrush(QColor("lightgray")));
  
  m_filesTable->setSelectionMode(QTableWidget::ExtendedSelection);
  string style = string("QTableWidget::indicator:unchecked ")
    + "{background-color:white; border: 1px solid black;}; " +
    "selection-background-color: rgba(128, 128, 128, 40);";

  m_filesTable->setSelectionMode(QTableWidget::NoSelection);

  m_filesTable->setStyleSheet(style.c_str());
  m_filesTable->resizeColumnsToContents();
  m_filesTable->resizeRowsToContents();
  
  //exec(); // Pop-up the filled in dialog

  // The processing of user's choice happens in aspView::showFilesChosenByUser()
  
  return;
}

