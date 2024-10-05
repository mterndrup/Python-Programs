

class Solution(object):

    #1. Two Sum
    def twoSum(self, nums, target):
        """
        :type nums: List[int]
        :type target: int
        :rtype: List[int]
        """

    #9. Palindrome Number
    def isPalindrome(self, x):
        """
        :type x: int
        :rtype: bool
        """
        sx = str(x)
        sx_reverse = sx[::-1]
        if (sx == sx_reverse):
            return True
        else:
            return False

    #13. Roman to Integer
    def romanToInt(self, s):
        """
        :type s: str
        :rtype: int
        """
        translations = {
            "I": 1,
            "V": 5,
            "X": 10,
            "L": 50,
            "C": 100,
            "D": 500,
            "M": 1000
        }
    
    #593. Valid Square
    def validSquare(self, p1, p2, p3, p4):
        """
        :type p1: List[int]
        :type p2: List[int]
        :type p3: List[int]
        :type p4: List[int]
        :rtype: bool
        """
        pass

sol = Solution()
print(sol.isPalindrome(int(input("Enter a number:"))))
        
