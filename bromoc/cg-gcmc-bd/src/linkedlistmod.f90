!--------------------------------------------------
!       Demonstration of linked lists
!
! Simple linked list of integers implemented to
! simulate a stack, i.e. it is operated in
! LIFO mode.
!
! (c) 1998 A. Migdalas, LiTH, Sweden
! (c) 2000 A. Migdalas @ MSCC
!
! 

MODULE Integer_Linked_List

  PRIVATE :: Node, Insert, Delete, SearchList, SearchListPosition

  TYPE Node
     INTEGER :: IntVal
     TYPE(Node), POINTER :: Next
  END TYPE Node

  TYPE IntList
     PRIVATE
     TYPE(Node), POINTER :: Header
  END TYPE IntList

  INTERFACE Push
     MODULE PROCEDURE Insert
  END INTERFACE

  INTERFACE Pop
     MODULE PROCEDURE Delete
  END INTERFACE

  INTERFACE Search
     MODULE PROCEDURE SearchList, SearchListPosition
  END INTERFACE

  CONTAINS

    !---------------------------------------------
    ! Creation of a new empty list.
    ! Usage:
    !       TYPE(IntList) :: listVar
    !       listVar = NewList()
    !
    
    FUNCTION NewList() RESULT(NewResult)
      IMPLICIT NONE

      TYPE(IntList) :: NewResult

      NULLIFY(NewResult%Header)

    END FUNCTION NewList

    !----------------------------------------------
    ! Reports empty list.
    ! Usage:
    !       TYPE(IntList) :: listVar
    !       IF ( EmptyList(listVar) ) THEN ...
    !

    LOGICAL FUNCTION EmptyList(listVar)
      IMPLICIT NONE

      TYPE(IntList), INTENT(IN) :: listVar
      
      EmptyList = .NOT.ASSOCIATED(listVar%Header)

    END FUNCTION EmptyList

    !------------------------------------------------
    ! Inserts a new item (here Integer) into the list.
    ! Insertion at the head of the list, i.e. operates
    ! like a stack.
    ! Usage:
    !       TYPE(IntList) :: listVar
    !       INTEGER       :: intItem
    !       CALL Insert(intItem,listVar)
    !

    SUBROUTINE Insert( intItem, listVar)
      IMPLICIT NONE

      TYPE(IntList), INTENT(INOUT) :: listVar
      INTEGER, INTENT(IN) :: intItem

      TYPE(Node), POINTER :: newNode

      ALLOCATE(newNode)
      newNode%IntVal = intItem

      IF ( EmptyList(listVar) ) THEN
         NULLIFY(newNode%Next)
         listVar%Header => newNode
      ELSE
         newNode%IntVal = intItem
         ALLOCATE(newNode%Next)
         newNode%Next => listVar%Header
         listVar%Header => newNode
      ENDIF

    END SUBROUTINE Insert

    !-----------------------------------------------
    ! Deletes an item (here Integer) from the list.
    ! Deletion from the head of the list, i.e. operates
    ! like a stack.
    ! Usage:
    !       TYPE(IntList) :: listVar
    !       INTEGER       :: intItem
    !       CALL Delete(intItem,listVar,notEmpty)
    !

    SUBROUTINE Delete( intItem, listVar, empty)
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: intItem
      TYPE(IntList), INTENT(INOUT) :: listVar
      LOGICAL, INTENT(OUT) :: empty

      TYPE(Node), POINTER :: helpNode

      empty = EmptyList(listVar)
      IF ( empty ) RETURN

      intItem = listVar%Header%IntVal

      ALLOCATE(helpNode)

      helpNode%Next => listVar%Header
      listVar%Header => helpNode%Next%Next

      DEALLOCATE(helpNode%Next, helpNode)

    END SUBROUTINE Delete

    !--------------------------------------------------
    ! Linearly Search the list for a specific value
    ! Usage:
    !       TYPE(IntList) :: listVar
    !       INTEGER       :: intItem
    !        IF ( SearchList(intItem,listVar) ) THEN ...
    !
    
    LOGICAL FUNCTION SearchList( intItem, listVar )
      IMPLICIT NONE

      TYPE(IntList), INTENT(IN) :: listVar
      INTEGER, INTENT(IN) :: intItem

      TYPE(Node), POINTER :: helpNode

      SearchList = EmptyList(listVar)

      IF ( SearchList ) RETURN

      helpNode => listVar%Header

      DO
         SearchList = ASSOCIATED(helpNode)
         IF ( .NOT.SearchList ) RETURN
         SearchList = ( intItem == helpNode%IntVal )
         IF ( SearchList ) RETURN
         helpNode => helpNode%Next
      ENDDO

    END FUNCTION SearchList

    !--------------------------------------------------
    ! Linearly Search the list for a specific value and
    ! return its "position" in the list if found
    ! Usage:
    !       TYPE(IntList) :: listVar
    !       INTEGER       :: intItem
    !       TYPE(Node)    :: position
    !       LOGICAL       :: found
    !       CALL SearchListPosition(intItem,listVar,position,found) 
    !       IF ( found ) THEN ...
    !
    
    SUBROUTINE SearchListPosition( intItem, listVar, posPointer, found )
      IMPLICIT NONE

      TYPE(IntList), INTENT(IN) :: listVar 
      TYPE(IntList), INTENT(OUT) :: posPointer
      INTEGER, INTENT(IN) :: intItem
      LOGICAL, INTENT(OUT) :: found             

      TYPE(Node), POINTER :: position

      found = .NOT.EmptyList(listVar)

      IF ( .NOT.found ) RETURN

      position => listVar%Header

      DO
         found = ASSOCIATED(position)
         IF ( .NOT.found ) RETURN
         found = ( intItem == position%IntVal )
         IF ( found ) THEN
            posPointer%Header => position
            RETURN
         ENDIF
         position => position%Next
      ENDDO

    END SUBROUTINE SearchListPosition

    !--------------------------------------------------
    ! Changes the value kept in a list-node pointed to by
    ! listPosition, which is a value returned by
    ! SearchListPosition.
    ! Usage:
    !       TYPE(IntList) :: listVar, listPosition
    !       INTEGER       :: intItem
    !       LOGICAL       :: found
    !       CALL SearchListPosition(intItem,listVar,listPosition,found) 
    !       CALL SetValue(intItem,listPosition)
    !

    SUBROUTINE SetValue( intItem, listPosition)
      IMPLICIT NONE

      TYPE(IntList), INTENT(INOUT) :: listPosition
      INTEGER, INTENT(IN) :: intItem

      listPosition%Header%IntVal = intItem

    END SUBROUTINE SetValue
         
      
      
    !---------------------------------------------------
    ! Prints the items in the list, if list is not empty.
    ! Usage:
    !       TYPE(IntList) :: listVar
    !       CALL PrintList(listVar)
    !
    
    SUBROUTINE PrintList(listVar)
      IMPLICIT NONE

      TYPE(IntList), INTENT(IN) :: listVar
      
      TYPE(Node), POINTER :: helpNode

      IF ( EmptyList(listVar) ) RETURN

      helpNode => listVar%Header

      DO
         IF( .NOT.ASSOCIATED(helpNode) ) EXIT
         PRINT *, helpNode%IntVal
         helpNode => helpNode%Next
      ENDDO

    END SUBROUTINE PrintList

         

  END MODULE Integer_Linked_List

