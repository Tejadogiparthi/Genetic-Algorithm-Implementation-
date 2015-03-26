Imports System.Math
Imports System
Imports System.Collections
Public Class Form1
    'All the variables that are declared are as follows
    Dim PopulationSize As Integer = 200 ' This specifies the population size here it is 200
    Dim Folding As Integer
    Public population(200) As genotype ' This is the first population
    Public newpopulation(200) As genotype 'This is the new population
    Dim ProteinStructure As String ' This stores the protein structure value that we give dynamically 
    Dim ProteinLength As Integer ' This stores the length of the above input input structure
    Dim HydrophobicPosition As Integer() 'This array is used to store all the indexes of hydropohobic positions from the input protein sequence
    Dim HydrophobicOccurences As Integer 'This stores the number of hydrophobic occurrences present in the input protein sequence
    Dim CurrentPosNewPopulation As Integer 'This keeps track of the index in the new population
    Dim CompleteFitness As Integer = 0
    Dim eliteRate As Decimal ' This stores the elite rate value that we give dynamically 
    Dim crossOverRate As Decimal ' This stores the cross over rate that we give dynamically 
    Dim mutationRate As Decimal ' This stores the mutation rate that we give dynamically
    Dim mutationPositionInNewPopulation As Integer
    Dim Generation As Integer = 0
    Dim maxGenerations As Integer ' This stores the maximum generations that we give dynamically 


    Public Class genotype
        Implements IComparable
        Public Fitness As Integer
        Public X(64) As Integer
        Public Y(64) As Integer

        'This function is used for sorting based on their fitness in ascending order
        Public Function CompareTo(ByVal gene As Object) As Integer _
        Implements IComparable.CompareTo
            If CType(gene, genotype).Fitness < Me.Fitness Then
                Return 1
            ElseIf CType(gene, genotype).Fitness = Me.Fitness Then
                Return 0
            ElseIf CType(gene, genotype).Fitness > Me.Fitness Then
                Return -1
            End If
            Return Nothing
        End Function
    End Class

    'This method checks the input protein structure and if it is hydrophoibc it stores it correspoing index into an array.
    Function Hyprophobicpositions()
        Dim HydrophobicIndex As Integer = 1
        ProteinLength = ProteinStructure.Length
        HydrophobicPosition = New Integer(ProteinLength) {}
        HydrophobicOccurences = 0
        Dim hOccurence As Char() = ProteinStructure.ToCharArray()
        For index = 1 To ProteinLength
            If (hOccurence(index - 1) = "h") Then
                HydrophobicPosition(HydrophobicIndex) = index
                HydrophobicIndex = HydrophobicIndex + 1
                HydrophobicOccurences = HydrophobicOccurences + 1
            End If
        Next index
        Return Nothing
    End Function

    Function Initialization()
        Dim i As Integer
        For i = 1 To PopulationSize
            Folding = 0
            RandomOrientation(i)

            While (Folding = 0)
                RandomOrientation(i)
            End While
            population(i).Fitness = ComputeFitness(i)
            CompleteFitness = CompleteFitness + population(i).Fitness
        Next i
        Return Nothing
    End Function

    

    'This will compute the fitness for a given protein structure
    Function ComputeFitness(n As Long) As Integer
        Dim isSequential As Integer
        Dim Fitness As Integer = 0
        Dim latticeDistance As Integer
        For i = 1 To HydrophobicOccurences - 1
            For j = i + 1 To HydrophobicOccurences
                isSequential = (Abs(HydrophobicPosition(i) - HydrophobicPosition(j))) '/*Not Sequential */
                If (isSequential <> 1) Then
                    latticeDistance = Abs(population(n).X(HydrophobicPosition(i)) - population(n).X(HydrophobicPosition(j))) + Abs(population(n).Y(HydrophobicPosition(i)) - population(n).Y(HydrophobicPosition(j)))
                    If (latticeDistance = 1) Then
                        Fitness = Fitness - 1
                    End If
                End If
            Next j
        Next i
        Return Fitness
    End Function

    ' This function will caluclate elite population
    Private Sub caluclatingElitePopulation()
        newpopulation = New genotype(PopulationSize) {}
        Dim elitePopulation As Integer
        elitePopulation = eliteRate * PopulationSize
        Array.ConstrainedCopy(population, 1, newpopulation, 1, elitePopulation)
    End Sub

    'This function will perform cross over to the population
    Private Sub CrossOverPopulation()
        Dim crossOverStartIndex As Integer = eliteRate * PopulationSize + 1
        Dim crossOverLastIndex As Integer = crossOverRate * PopulationSize + crossOverStartIndex - 1
        Dim crossOverPoint As Integer
        Dim i, j As Integer
        Dim maxEndPoint As Integer = ProteinLength - 3
        For Index = crossOverStartIndex To crossOverLastIndex
            CurrentPosNewPopulation = Index
            Do Until i > 0
                i = ChromosomeSelectionUsingRoulettewheelSelection()
            Loop
            Do Until j > 0
                j = ChromosomeSelectionUsingRoulettewheelSelection()
            Loop
            newpopulation(CurrentPosNewPopulation) = New genotype()
            Randomize()
            crossOverPoint = (maxEndPoint * Rnd() + 2)
            Dim Success As Integer = CrossOver(i, j, crossOverPoint)
            While Success = 0
                Do Until i > 0
                    i = ChromosomeSelectionUsingRoulettewheelSelection()
                Loop
                Do Until j > 0
                    j = ChromosomeSelectionUsingRoulettewheelSelection()
                Loop
                Randomize()
                crossOverPoint = (maxEndPoint * Rnd() + 2)
                Success = CrossOver(i, j, crossOverPoint)
            End While
        Next Index
    End Sub

    'This function will fill the remaining population in new population array
    Private Sub FillRemainingNewPopulation()
        Try
            Dim remainingNewPopulationStartIndex As Integer = eliteRate * PopulationSize + crossOverRate * PopulationSize + 1
            Array.ConstrainedCopy(population, remainingNewPopulationStartIndex, newpopulation, remainingNewPopulationStartIndex, PopulationSize - remainingNewPopulationStartIndex + 1)
        Catch ex As Exception
            MessageBox.Show(ex.Message)
        End Try
    End Sub

    Function RandomOrientation(m As Long)

        Dim PreviousDirection, PresentDirection, i, temp1, temp2, temp3, X, Y, j, Flag, Step2 As Integer
        Dim a(4), Ax(4), Ay(4) As Integer

        '                                        3
        '             Select Direction as:     2 X 1
        '                                        4
        '
        population(m) = New genotype()
        Folding = 1
        population(m).X(1) = 0
        population(m).Y(1) = 0
        population(m).X(2) = 1
        population(m).Y(2) = 0
        PreviousDirection = 1


        For i = 3 To ProteinLength

            Select Case PreviousDirection
                Case 1
                    a(1) = 1
                    Ax(1) = 1
                    Ay(1) = 0
                    a(2) = 3
                    Ax(2) = 0
                    Ay(2) = 1
                    a(3) = 4
                    Ax(3) = 0
                    Ay(3) = -1
                Case 2
                    a(1) = 2
                    Ax(1) = -1
                    Ay(1) = 0
                    a(2) = 3
                    Ax(2) = 0
                    Ay(2) = 1
                    a(3) = 4
                    Ax(3) = 0
                    Ay(3) = -1
                Case 3
                    a(1) = 1
                    Ax(1) = 1
                    Ay(1) = 0
                    a(2) = 2
                    Ax(2) = -1
                    Ay(2) = 0
                    a(3) = 3
                    Ax(3) = 0
                    Ay(3) = 1
                Case 4
                    a(1) = 1
                    Ax(1) = 1
                    Ay(1) = 0
                    a(2) = 2
                    Ax(2) = -1
                    Ay(2) = 0
                    a(3) = 4
                    Ax(3) = 0
                    Ay(3) = -1
            End Select

            temp1 = Int(3 * Rnd() + 1)
            PresentDirection = temp1
            temp2 = 0
            temp3 = 0
            X = population(m).X(i - 1) + Ax(temp1)
            Y = population(m).Y(i - 1) + Ay(temp1)
            Flag = 0

            For j = 1 To i - 1
                If (X = population(m).X(j) And Y = population(m).Y(j)) Then
                    Flag = 1
                    GoTo MyJump1
                End If
            Next j

MyJump1:
            If (Flag = 1) Then
                Flag = 0
                Step2 = 6 - temp1
                Select Case Step2
                    Case 3
                        If Int(Rnd() * 2 + 1) = 1 Then
                            temp2 = 1
                        Else
                            temp2 = 2
                        End If
                    Case 4
                        If Int(Rnd() * 2 + 1) = 1 Then
                            temp2 = 1
                        Else
                            temp2 = 3
                        End If
                    Case 5
                        If Int(Rnd() * 2 + 1) = 1 Then
                            temp2 = 2
                        Else
                            temp2 = 3
                        End If
                End Select

                PresentDirection = temp2
                temp3 = 6 - (temp1 + temp2)
                X = population(m).X(i - 1) + Ax(temp2)
                Y = population(m).Y(i - 1) + Ay(temp2)

                For j = 1 To i - 1
                    If (X = population(m).X(j) And Y = population(m).Y(j)) Then
                        Flag = 1
                        GoTo MyJump2
                    End If
                Next j
MyJump2:
                If (Flag = 1) Then
                    Flag = 0
                    PresentDirection = temp3
                    X = population(m).X(i - 1) + Ax(temp3)
                    Y = population(m).Y(i - 1) + Ay(temp3)
                    For j = 1 To i - 1
                        If (X = population(m).X(j) And Y = population(m).Y(j)) Then
                            Flag = 1
                            Folding = 0
                            'GoTo MyJump3

                        End If
                    Next j
                End If
            End If
            PreviousDirection = a(PresentDirection)
            population(m).X(i) = population(m).X(i - 1) + Ax(PresentDirection)
            population(m).Y(i) = population(m).Y(i - 1) + Ay(PresentDirection)
        Next i
MyJump3:
        Return Nothing
    End Function

    'This method selects chromosome for crossover using  the roulettewheel selection logic
    Private Function ChromosomeSelectionUsingRoulettewheelSelection() As Integer
        Randomize()
        Dim rndVar As Integer = Rnd() * Math.Abs(CompleteFitness)
        Dim index As Integer
        For index = 1 To PopulationSize
            rndVar = rndVar - Math.Abs(population(index).Fitness)
            If (rndVar < 0) Then
                Return index - 1
            End If
        Next
        Return Nothing
    End Function
    Function CrossOver(i As Long, j As Long, n As Integer) As Long

        Dim PrevDirection, k, z, p As Long
        Dim temp1, temp2, temp3, Collision, dx, dy, Step2 As Long
        Dim id As Long
        Dim a(0 To 4) As Long
        Dim Ax(0 To 4) As Long
        Dim Ay(0 To 4) As Long

        id = CurrentPosNewPopulation

        '/* Detect Previous Direction */
        If (population(i).X(n) = population(i).X(n - 1)) Then
            p = population(i).Y(n - 1) - population(i).Y(n)
            If (p = 1) Then
                PrevDirection = 3
            Else
                PrevDirection = 4
            End If

        Else
            p = population(i).X(n - 1) - population(i).X(n)
            If (p = 1) Then
                PrevDirection = 1
            Else
                PrevDirection = 2
            End If
        End If


        Select Case PrevDirection
            Case 1
                Ax(1) = -1
                Ay(1) = 0
                Ax(2) = 0
                Ay(2) = 1
                Ax(3) = 0
                Ay(3) = -1
            Case 2
                Ax(1) = 1
                Ay(1) = 0
                Ax(2) = 0
                Ay(2) = 1
                Ax(3) = 0
                Ay(3) = -1
            Case 3
                Ax(1) = 1
                Ay(1) = 0
                Ax(2) = -1
                Ay(2) = 0
                Ax(3) = 0
                Ay(3) = -1

            Case 4
                Ax(1) = 1
                Ay(1) = 0
                Ax(2) = -1
                Ay(2) = 0
                Ax(3) = 0
                Ay(3) = 1
        End Select

        temp1 = Int(Rnd() * 3 + 1)

        newpopulation(id).X(n + 1) = population(i).X(n) + Ax(temp1)
        newpopulation(id).Y(n + 1) = population(i).Y(n) + Ay(temp1)
        Collision = 0

        dx = newpopulation(id).X(n + 1) - population(j).X(n + 1)
        dy = newpopulation(id).Y(n + 1) - population(j).Y(n + 1)

        For k = n + 1 To ProteinLength
            newpopulation(id).X(k) = population(j).X(k) + dx

            newpopulation(id).Y(k) = population(j).Y(k) + dy

            For z = 1 To n
                If ((newpopulation(id).X(k) = population(i).X(z)) And (newpopulation(id).Y(k) = population(i).Y(z))) Then
                    Collision = 1
                    GoTo MyOut1
                End If
            Next z
        Next k

MyOut1:
        If (Collision = 1) Then         '/* ======> Second try ==== */
            Collision = 0
            Step2 = 6 - temp1
            Select Case Step2
                Case 3
                    If Int(Rnd() * 2 + 1) = 1 Then
                        temp2 = 1
                    Else
                        temp2 = 2
                    End If

                Case 4
                    If Int(Rnd() * 2 + 1) = 1 Then
                        temp2 = 1
                    Else
                        temp2 = 3
                    End If

                Case 5
                    If Int(Rnd() * 2 + 1) = 1 Then
                        temp2 = 2
                    Else
                        temp2 = 3
                    End If
            End Select

            temp3 = 6 - (temp1 + temp2)
            newpopulation(id).X(n + 1) = population(i).X(n) + Ax(temp2)
            newpopulation(id).Y(n + 1) = population(i).Y(n) + Ay(temp2)
            dx = newpopulation(id).X(n + 1) - population(j).X(n + 1)
            dy = newpopulation(id).Y(n + 1) - population(j).Y(n + 1)

            For k = n + 1 To ProteinLength

                newpopulation(id).X(k) = population(j).X(k) + dx
                newpopulation(id).Y(k) = population(j).Y(k) + dy

                For z = 1 To n
                    If ((newpopulation(id).X(k) = population(i).X(z)) And (newpopulation(id).Y(k) = population(i).Y(z))) Then
                        Collision = 1
                        GoTo MyOut2
                    End If
                Next z
            Next k

MyOut2:
            If (Collision = 1) Then
                Collision = 0
                newpopulation(id).X(n + 1) = population(i).X(n) + Ax(temp3)
                newpopulation(id).Y(n + 1) = population(i).Y(n) + Ay(temp3)
                dx = newpopulation(id).X(n + 1) - population(j).X(n + 1)
                dy = newpopulation(id).Y(n + 1) - population(j).Y(n + 1)
                For k = n + 1 To ProteinLength
                    newpopulation(id).X(k) = population(j).X(k) + dx
                    newpopulation(id).Y(k) = population(j).Y(k) + dy
                    For z = 1 To n
                        If ((newpopulation(id).X(k) = population(i).X(z)) And (newpopulation(id).Y(k) = population(i).Y(z))) Then
                            Collision = 1
                            GoTo MyOut3
                        End If
                    Next z
                Next k
            End If '/* 3rd try if ends */
        End If '/* 2nd try if ends */

MyOut3:
        If Collision = 0 Then
            For k = 1 To n
                newpopulation(id).X(k) = population(i).X(k)
                newpopulation(id).Y(k) = population(i).Y(k)
            Next k
            CrossOver = 1
        End If

    End Function

    'This function will perform mutation
    Private Sub PerformMutation()

        Dim mutationPopulation As Integer = mutationRate * PopulationSize
        Randomize()
        Dim geneToBeMutated As Integer = 199 * Rnd() + 1
        Randomize()
        Dim maxEndPoint As Integer = ProteinLength - 3
        Dim mutationPoint As Integer = maxEndPoint * Rnd() + 2
        Try
            Randomize()
            mutationPositionInNewPopulation = 189 * Rnd() + 11
            For index = 1 To mutationPopulation
                mutationPositionInNewPopulation = mutationPositionInNewPopulation
                Dim MutationStatus As Integer = Mutation(geneToBeMutated, mutationPoint)
                While MutationStatus = 0
                    geneToBeMutated = 199 * Rnd() + 1
                    mutationPoint = maxEndPoint * Rnd() + 2
                    MutationStatus = Mutation(geneToBeMutated, mutationPoint)
                End While
            Next
        Catch ex As Exception
            MessageBox.Show(geneToBeMutated + "geneToBeMutated" + mutationPoint + "mutationPoint" + mutationPositionInNewPopulation)
        End Try
    End Sub

    Function Mutation(i As Long, n As Integer) As Long
        Dim id As Long
        Dim a As Long
        Dim b As Long
        Dim A_Limit As Long
        Dim choice As Long
        Dim Collision As Long
        Dim k As Long
        Dim z As Long
        Dim p As Long
        Dim Ary(3) As Integer

        id = mutationPositionInNewPopulation

        ' possible rotations 90ß,180ß,270ß
        '           index       1   2    3
        '


        Ary(1) = 1
        Ary(2) = 2
        Ary(3) = 3
        A_Limit = 3

        a = population(i).X(n)          '/* (a, b) rotating point */
        b = population(i).Y(n)

        Do
            Collision = 0
            If (A_Limit > 1) Then
                Randomize()
                choice = Int(A_Limit * Rnd() + 1)
            Else
                choice = A_Limit
            End If


            p = Ary(choice)
            For k = choice To A_Limit - 1
                Ary(k) = Ary(k + 1)
            Next k

            A_Limit = A_Limit - 1

            For k = n + 1 To ProteinLength
                Select Case p

                    Case 1
                        newpopulation(id).X(k) = a + b - population(i).Y(k)       '/* X' = (a+b)-Y  */
                        newpopulation(id).Y(k) = population(i).X(k) + b - a       '/* Y' = (X+b)-a  */
                    Case 2
                        newpopulation(id).X(k) = 2 * a - population(i).X(k)       '/* X' = (2a - X) */
                        newpopulation(id).Y(k) = 2 * b - population(i).Y(k)       '/* Y' = (2b-Y)   */
                    Case 3
                        newpopulation(id).X(k) = population(i).Y(k) + a - b       '/* X' =  Y+a-b   */
                        newpopulation(id).Y(k) = a + b - population(i).X(k)       '/* Y' =  (a+b)-X */
                End Select

                For z = 1 To n

                    If ((newpopulation(id).X(k) = population(i).X(z)) And (newpopulation(id).Y(k) = population(i).Y(z))) Then
                        Collision = 1
                        GoTo MyJump
                    End If
                Next z
            Next k

            If (Collision = 0) Then
                A_Limit = 0
            End If
MyJump:
        Loop Until A_Limit = 0

        If (Collision = 0) Then
            For k = 1 To n
                newpopulation(id).X(k) = population(i).X(k)
                newpopulation(id).Y(k) = population(i).Y(k)
            Next k
            Mutation = 1
        Else
            Mutation = 0
        End If

    End Function

    'This function will compute the next generation
    Private Sub NextGeneration()
        Do
            Chart1.Series("Protein Structure").Points.Clear()
            CompleteFitness = 0
            Array.ConstrainedCopy(newpopulation, 1, population, 1, PopulationSize)
            For index = 1 To PopulationSize
                population(index).Fitness = 0
                population(index).Fitness = ComputeFitness(index)
                CompleteFitness = CompleteFitness + population(index).Fitness
            Next

            'Sorting the Population based on the Fitness 
            Array.Sort(population)

            'Caluclate the elite population
            caluclatingElitePopulation()

            'CrossOver Population
            CrossOverPopulation()

            'Fill Remaining Population
            FillRemainingNewPopulation()

            ' perform Mutation
            PerformMutation()
            Generation = Generation + 1
            TextBox6.Text = population(1).Fitness
            TextBox6.Enabled = False
            Chart1.Series("Protein Structure").BorderWidth = 5
            Chart1.ChartAreas(0).AxisX.Interval = 1
            Chart1.ChartAreas(0).AxisY.Interval = 1
            For index = 1 To ProteinStructure.Length
                Chart1.Series("Protein Structure").Points.AddXY(population(1).X(index), population(1).Y(index))
            Next
            For c = 0 To ProteinLength - 1
                'HydropPhobic Positions is represented in Black Color
                If HydrophobicPosition.Contains(c + 1) Then
                    Chart1.Series("Protein Structure").Points(c).MarkerStyle = DataVisualization.Charting.MarkerStyle.Circle
                    Chart1.Series("Protein Structure").Points(c).MarkerSize = 10
                    Chart1.Series("Protein Structure").Points(c).MarkerColor = Color.Black
                Else
                    'HydropPhobic Positions is represented in Red Color
                    Chart1.Series("Protein Structure").Points(c).MarkerStyle = DataVisualization.Charting.MarkerStyle.Circle
                    Chart1.Series("Protein Structure").Points(c).MarkerSize = 10
                    Chart1.Series("Protein Structure").Points(c).MarkerColor = Color.Red
                End If
            Next
            TextBox7.Text = Generation
            TextBox7.Enabled = False
            Chart1.Refresh()
            TextBox6.Refresh()
            TextBox7.Refresh()
            System.Threading.Thread.Sleep(30)
        Loop Until Generation = maxGenerations
        Return
    End Sub

    'This function will start  to execute after taking the inputs dynamically
    Private Sub StartExecution()
        Hyprophobicpositions()
        Initialization()
        Array.Sort(population)
        caluclatingElitePopulation()
        CrossOverPopulation()
        FillRemainingNewPopulation()
        PerformMutation()
        NextGeneration()
    End Sub

    ' This function will get executed when start button clicked after taking the inputs dynamically
    Private Sub Button1_Click(sender As Object, e As EventArgs) Handles Button1.Click
        Try
            Generation = 0
            ProteinStructure = TextBox1.Text
            eliteRate = TextBox2.Text
            crossOverRate = TextBox3.Text
            mutationRate = TextBox4.Text
            maxGenerations = Val(TextBox5.Text)
            StartExecution()
        Catch ex As Exception
            MessageBox.Show(ex.Message)
        End Try

    End Sub

    ' This function will get executed when reset button is clicked 
    Private Sub Button2_Click(sender As Object, e As EventArgs) Handles Button2.Click
        TextBox1.ResetText()
        TextBox2.ResetText()
        TextBox3.ResetText()
        TextBox4.ResetText()
        TextBox5.ResetText()
        TextBox6.ResetText()
        TextBox7.ResetText()
        Chart1.Series("Protein Structure").Points.Clear()
        Generation = 0
    End Sub




End Class
